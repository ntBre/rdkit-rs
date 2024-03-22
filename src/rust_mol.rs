//! Experimental module for a Rust version of [ROMol][crate::ROMol] constructed
//! using serde and [ROMol::to_json][crate::ROMol::to_json].

use serde::{Deserialize, Serialize};

use crate::{RDError, ROMol};

#[derive(Serialize, Deserialize)]
struct RDKitJSON {
    version: usize,
}

#[derive(Serialize, Deserialize)]
struct AtomDefaults {
    z: usize,
    #[serde(rename = "impHs")]
    imp_hs: usize,
    chg: usize,
    #[serde(rename = "nRad")]
    n_rad: usize,
    isotope: usize,
    stereo: String,
}

#[derive(Serialize, Deserialize)]
struct BondDefaults {
    bo: usize,
    stereo: String,
}

#[derive(Serialize, Deserialize)]
struct Defaults {
    atom: AtomDefaults,
    bond: BondDefaults,
}

#[derive(Serialize, Deserialize)]
struct Atom {
    #[serde(skip_serializing_if = "Option::is_none")]
    z: Option<usize>,
    #[serde(rename = "impHs", skip_serializing_if = "Option::is_none")]
    imp_hs: Option<usize>,
}

#[derive(Serialize, Deserialize)]
struct Bond {
    #[serde(skip_serializing_if = "Option::is_none")]
    bo: Option<usize>,
    atoms: Vec<usize>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
struct Extension {
    name: String,
    format_version: usize,
    toolkit_version: String,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    aromatic_atoms: Vec<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    aromatic_bonds: Vec<usize>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    atom_rings: Vec<Vec<usize>>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    cip_ranks: Vec<usize>,
}

#[derive(Serialize, Deserialize)]
struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    extensions: Vec<Extension>,
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RSMol {
    rdkitjson: RDKitJSON,
    defaults: Defaults,
    molecules: Vec<Molecule>,
}

impl RSMol {
    pub fn from_json(s: &str) -> Result<Self, RDError> {
        Ok(serde_json::from_str(s)?)
    }

    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap()
    }
}

impl From<ROMol> for RSMol {
    fn from(value: ROMol) -> Self {
        Self::from_json(&value.to_json()).unwrap()
    }
}

impl From<RSMol> for ROMol {
    fn from(value: RSMol) -> Self {
        Self::from_json(&value.to_json())
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::read_to_string, path::PathBuf};

    use super::*;

    #[test]
    fn test_from_json() {
        let dir = PathBuf::from("testfiles");
        for file in ["rdkit", "romol"] {
            eprintln!("testing {file}.json");
            let f = dir.join(file).with_extension("json");
            let s = read_to_string(f).unwrap();
            let got = RSMol::from_json(&s).unwrap().to_json();
            assert_eq!(got, s.trim(), "got = {got}, want = {}", s.trim());
        }
    }

    #[test]
    fn test_from_romol() {
        let mol = ROMol::from_smiles("CCO");
        let _ = RSMol::from(mol);
    }

    #[test]
    fn test_into_romol() {
        let s = read_to_string("testfiles/rdkit.json").unwrap();
        let mol = RSMol::from_json(&s).unwrap();
        let _ = ROMol::from(mol);
    }
}
