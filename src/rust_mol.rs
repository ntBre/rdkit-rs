//! Experimental module for a Rust version of [ROMol][crate::ROMol] constructed
//! using serde and [ROMol::to_json][crate::ROMol::to_json].

use serde::{Deserialize, Serialize};

use crate::RDError;

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
    aromatic_atoms: Vec<usize>,
    aromatic_bonds: Vec<usize>,
    atom_rings: Vec<Vec<usize>>,
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

#[cfg(test)]
mod tests {
    use std::fs::read_to_string;

    use super::*;

    #[test]
    fn test_from_json() {
        let s = read_to_string("testfiles/rdkit.json").unwrap();
        let got = RSMol::from_json(&s).unwrap().to_json();
        assert_eq!(got, s.trim());
    }
}
