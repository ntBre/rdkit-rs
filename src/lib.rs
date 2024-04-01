use std::{
    collections::HashMap,
    ffi::{c_int, c_uint, CStr, CString},
    fmt::Display,
};

use bitflags::bitflags;
use rdkit_sys::{
    RDKit_JSONToMol, RDKit_MolToJSON, RDKit_MolToSmiles, RDKit_ROMol,
    RDKit_ROMol_delete, RDKit_SmartsToMol, RDKit_SmilesToMol,
};

use self::bitvector::BitVector;

#[cfg(test)]
mod tests;

pub mod bitvector;
pub mod errors;
pub mod fingerprint;
pub mod fragment;
pub mod mol_supplier;
pub mod rust_mol;

pub use errors::RDError;
pub use mol_supplier::SDMolSupplier;
pub use rdkit_sys::Point3D;

pub struct SmilesParserParams {
    /// defaults to true
    pub remove_hs: bool,
    /// defaults to true
    pub sanitize: bool,
}

impl Default for SmilesParserParams {
    fn default() -> Self {
        Self {
            remove_hs: true,
            sanitize: true,
        }
    }
}

pub struct ROMol(*mut RDKit_ROMol);

pub struct Conformer(*mut rdkit_sys::RDKit_Conformer);

impl Drop for Conformer {
    fn drop(&mut self) {
        unsafe {
            rdkit_sys::RDKit_Conformer_delete(self.0);
        }
    }
}

impl Conformer {
    pub fn get_positions(&self) -> Vec<Point3D> {
        unsafe {
            let mut n = 0;
            let ps = rdkit_sys::RDKit_Conformer_getPositions(self.0, &mut n);
            assert!(!ps.is_null());
            Vec::from_raw_parts(ps, n as usize, n as usize)
        }
    }
}

impl Display for ROMol {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_smiles())
    }
}

impl ROMol {
    pub fn from_smiles(smiles: &str) -> Self {
        Self::from_smiles_full(smiles, SmilesParserParams::default())
    }

    pub fn from_smiles_full(smiles: &str, params: SmilesParserParams) -> Self {
        let s = CString::new(smiles).expect("failed to create CString");
        unsafe {
            Self(RDKit_SmilesToMol(
                s.as_ptr(),
                params.remove_hs,
                params.sanitize,
            ))
        }
    }

    pub fn from_smarts(smarts: &str) -> Self {
        let s = CString::new(smarts).expect("failed to create CString");
        unsafe { Self(RDKit_SmartsToMol(s.as_ptr())) }
    }

    /// Create an [ROMol] from a JSON string. The format can be either
    /// [CommonChem](https://github.com/CommonChem/CommonChem), or the RDKit
    /// [extension](http://rdkit.org/docs/source/rdkit.Chem.rdMolInterchange.html).
    pub fn from_json(json: &str) -> Self {
        let s = CString::new(json).expect("failed to create CString");
        unsafe { Self(RDKit_JSONToMol(s.as_ptr())) }
    }

    /// Convert the molecule to RDKit's extension to the commonchem JSON format.
    pub fn to_json(&self) -> String {
        unsafe {
            let json = RDKit_MolToJSON(self.0);
            let s = CStr::from_ptr(json);
            s.to_str().unwrap().to_owned()
        }
    }

    pub fn to_smiles(&self) -> String {
        unsafe {
            let smiles = RDKit_MolToSmiles(self.0);
            let s = CStr::from_ptr(smiles);
            s.to_str().unwrap().to_owned()
        }
    }

    pub fn to_inchi_key(&self) -> String {
        unsafe {
            let smiles = rdkit_sys::RDKit_MolToInchiKey(self.0);
            let s = CStr::from_ptr(smiles);
            s.to_str().unwrap().to_owned()
        }
    }

    pub fn num_atoms(&self) -> usize {
        unsafe { rdkit_sys::RDKit_ROMol_getNumAtoms(self.0) as usize }
    }

    pub fn elements(&self) -> Vec<usize> {
        unsafe {
            let mut natoms = 0;
            let ptr = rdkit_sys::RDKit_ROMol_getElements(self.0, &mut natoms);
            let ret = Vec::from_raw_parts(ptr, natoms, natoms);
            ret.into_iter().map(|i| i as usize).collect()
        }
    }

    pub fn sanitize(&mut self, ops: SanitizeFlags) {
        unsafe {
            let status = rdkit_sys::RDKit_SanitizeMol(self.0, ops.bits());
            if status != 0 {
                panic!("sanitization failed");
            }
        }
    }

    pub fn set_aromaticity(&mut self, mdl: AromaticityModel) {
        unsafe {
            rdkit_sys::RDKit_SetAromaticity(self.0, mdl.bits());
        }
    }

    pub fn assign_stereochemistry(&mut self) {
        unsafe {
            rdkit_sys::RDKit_AssignStereochemistry(self.0);
        }
    }

    pub fn add_hs(&mut self) {
        unsafe {
            rdkit_sys::RDKit_AddHs(self.0);
        }
    }

    /// performs the sequence of cleaning operations used by the openff-toolkit:
    /// sanitize the molecule with ALL ^ ADJUSTHS ^ SETAROMATICITY, set
    /// aromaticity to MDL, assign stereochemistry, and add hydrogens.
    pub fn openff_clean(&mut self) {
        self.sanitize(
            SanitizeFlags::ALL
                ^ SanitizeFlags::ADJUSTHS
                ^ SanitizeFlags::SETAROMATICITY,
        );
        self.set_aromaticity(AromaticityModel::MDL);
        self.assign_stereochemistry();
        self.add_hs();
    }

    pub fn morgan_fingerprint(&self, radius: c_uint) -> HashMap<usize, usize> {
        unsafe {
            let mut len = 0;
            let ptr =
                rdkit_sys::RDKit_MorganFingerprint(self.0, radius, &mut len);
            let v1 = Vec::from_raw_parts(ptr, len, len);
            v1.into_iter()
                .map(|v| (v.bit as usize, v.count as usize))
                .collect()
        }
    }

    pub fn morgan_fingerprint_bit_vec<const N: usize>(
        &self,
        radius: c_uint,
    ) -> BitVector {
        unsafe {
            let mut tmp = [false; N];
            rdkit_sys::RDKit_MorganFingerprintBitVector(
                self.0,
                radius,
                N,
                tmp.as_mut_ptr(),
            );
            tmp.as_ref().into()
        }
    }

    pub fn draw_svg(
        &self,
        width: usize,
        height: usize,
        legend: &str,
        highlight_atoms: &[usize],
    ) -> String {
        let atoms: Vec<_> =
            highlight_atoms.iter().map(|a| *a as c_int).collect();
        unsafe {
            let legend = CString::new(legend).unwrap();
            let s = rdkit_sys::RDKit_MolDrawSVG(
                self.0,
                width as c_int,
                height as c_int,
                legend.as_ptr(),
                atoms.as_ptr(),
                highlight_atoms.len(),
            );
            CString::from_raw(s).to_str().unwrap().to_owned()
        }
    }

    pub fn replace_substructs(
        &self,
        query: &ROMol,
        replacement: &ROMol,
        all: bool,
    ) -> Vec<ROMol> {
        unsafe {
            let mut len = 0;
            let res = rdkit_sys::RDKit_ReplaceSubstructs(
                self.0,
                query.0,
                replacement.0,
                all,
                &mut len,
            );
            assert!(!res.is_null());
            let mut ret = Vec::with_capacity(len);
            for i in 0..len {
                ret.push(ROMol(*res.add(i)));
            }
            ret
        }
    }

    pub fn compute_2d_coords(&self) -> usize {
        unsafe { rdkit_sys::RDKit_compute2DCoords(self.0) as usize }
    }

    pub fn get_conformer(&self, id: usize) -> Conformer {
        unsafe {
            Conformer(rdkit_sys::RDKit_ROMol_getConformer(self.0, id as c_int))
        }
    }

    /// Call the sequence of operations to generate and return a vector of 2D
    /// coordinates
    pub fn get_2d_coords(&self) -> Vec<Point3D> {
        self.get_conformer(self.compute_2d_coords()).get_positions()
    }
}

impl Clone for ROMol {
    fn clone(&self) -> Self {
        unsafe { Self(rdkit_sys::RDKit_ROMol_copy(self.0)) }
    }
}

impl Drop for ROMol {
    fn drop(&mut self) {
        unsafe {
            RDKit_ROMol_delete(self.0);
        }
    }
}

unsafe impl Send for ROMol {}
unsafe impl Sync for ROMol {}

bitflags! {
    pub struct SanitizeFlags: c_uint {
        const NONE =                    0x0;
        const CLEANUP =                 0x1;
        const PROPERTIES =              0x2;
        const SYMMRINGS =               0x4;
        const KEKULIZE =                0x8;
        const FINDRADICALS =            0x10;
        const SETAROMATICITY =          0x20;
        const SETCONJUGATION =          0x40;
        const SETHYBRIDIZATION =        0x80;
        const CLEANUPCHIRALITY =        0x100;
        const ADJUSTHS =                0x200;
        const CLEANUP_ORGANOMETALLICS = 0x400;
        const ALL =                     0xFFFFFFF;
    }
}

bitflags! {
    pub struct AromaticityModel: c_uint {
        const DEFAULT = 0x0;
        const RDKIT = 0x1;
        const SIMPLE = 0x2;
        const MDL = 0x4;
        const CUSTOM = 0xFFFFFFF;
    }
}

pub fn find_smarts_matches(mol: &ROMol, smarts: &str) -> Vec<Vec<usize>> {
    let mut len = 0;
    let mut match_size = 0;
    let smarts = CString::new(smarts).unwrap();
    unsafe {
        let matches = rdkit_sys::find_smarts_matches(
            mol.0,
            smarts.as_ptr(),
            &mut len,
            &mut match_size,
        );
        let matches = Vec::from_raw_parts(matches, len, len);

        let mut ret = Vec::new();
        for mat in matches.chunks(match_size) {
            ret.push(mat.iter().map(|&x| x as usize).collect());
        }
        ret
    }
}

/// returns the sequence of chemical environment "tuples" that match `smarts` in
/// `mol`.
pub fn find_smarts_matches_mol(mol: &ROMol, smarts: &ROMol) -> Vec<Vec<usize>> {
    let mut len = 0;
    let mut match_size = 0;
    unsafe {
        let matches = rdkit_sys::find_smarts_matches_mol(
            mol.0,
            smarts.0,
            &mut len,
            &mut match_size,
        );
        let matches = Vec::from_raw_parts(matches, len, len);

        let mut ret = Vec::new();
        for mat in matches.chunks(match_size) {
            ret.push(mat.iter().map(|&x| x as usize).collect());
        }
        ret
    }
}
