use std::ffi::CString;

use rdkit_sys::{
    RDKit_SDMolSupplier, RDKit_create_mol_supplier, RDKit_delete_mol_supplier,
    RDKit_mol_supplier_at_end, RDKit_mol_supplier_next,
};

use crate::{RDError, ROMol};

pub mod multithreaded;

pub struct SDMolSupplier(*mut RDKit_SDMolSupplier);

unsafe impl Send for SDMolSupplier {}

impl SDMolSupplier {
    /// construct an [SDMolSupplier] from a filepath that can be converted to a
    /// CString. panics if this conversion fails.
    pub fn new(path: impl Into<Vec<u8>>) -> Result<Self, RDError> {
        let cpath = CString::new(path).expect("failed to create CString");
        unsafe {
            let inner = RDKit_create_mol_supplier(cpath.as_ptr(), true);
            if inner.is_null() {
                return Err(RDError);
            }
            Ok(Self(inner))
        }
    }

    /// reports whether or not `self` is at the end of the underlying file
    pub fn at_end(&self) -> bool {
        unsafe { RDKit_mol_supplier_at_end(self.0) }
    }
}

impl Drop for SDMolSupplier {
    fn drop(&mut self) {
        unsafe {
            RDKit_delete_mol_supplier(self.0);
        }
    }
}

impl Iterator for SDMolSupplier {
    type Item = Result<ROMol, RDError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.at_end() {
            return None;
        }
        unsafe {
            let mol = RDKit_mol_supplier_next(self.0);
            if mol.is_null() {
                return Some(Err(RDError));
            }
            Some(Ok(ROMol(mol)))
        }
    }
}
