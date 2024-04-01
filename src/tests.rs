use std::fs::read_to_string;

use super::*;

#[test]
fn to_inchi_key() {
    let benzene = ROMol::from_smiles("C1=CC=CC=C1");
    let got = benzene.to_inchi_key();
    let want = "UHOVQNZJYSORNB-UHFFFAOYSA-N";
    assert_eq!(got, want);
}

#[test]
fn rdkit_json() {
    let s = read_to_string("testfiles/rdkit.json").unwrap();
    ROMol::from_json(&s);
}

#[test]
fn commonchem_json() {
    let s = read_to_string("testfiles/commonchem.json").unwrap();
    ROMol::from_json(&s);
}

#[test]
fn elements() {
    let mol = ROMol::from_smiles("CCO");
    let got = mol.elements();
    let want = [6, 6, 8];
    assert_eq!(got, want);
}

#[test]
fn get_2d_coords() {
    let mol = ROMol::from_smiles("CCO");
    let coords = mol.get_2d_coords();
    assert_eq!(coords.len(), 3);
}
