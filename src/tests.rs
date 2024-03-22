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
