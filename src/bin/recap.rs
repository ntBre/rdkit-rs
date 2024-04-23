use rdkit_rs::{fragment::recap_decompose, ROMol};

// the smiles in main right now took ~100 ms in my chembl code. this one
// took 50.4 seconds

// [H]OC([H])([H])[C@]1([H])O[C@]([H])(O[C@@]2([H])[C@]([H])(C([H])([H])[H])O[C@@]([H])(C([H])([H])OC([H])([H])[C@@]3([H])[C@@]([H])(OC([H])([H])[H])C([H])([H])[C@]([H])(C([H])([H])OC([H])([H])[C@]4([H])[C@]([H])(C([H])([H])[H])O[C@@]([H])(C([H])([H])OC([H])([H])[C@@]5([H])[C@@]([H])(OC([H])([H])[H])C([H])([H])[C@]([H])(O[C@]6([H])C([H])([H])C7=C([H])C([H])([H])[C@]8(O[H])[C@]([H])(C([H])([H])[C@@]([H])(OC(=O)/C([H])=C(\[H])c9c([H])c([H])c([H])c([H])c9[H])[C@@]9(C([H])([H])[H])[C@]8(O[H])C([H])([H])C([H])([H])[C@@]9(O[H])C(=O)C([H])([H])[H])[C@@]7(C([H])([H])[H])C([H])([H])C6([H])[H])O[C@]5([H])C([H])([H])[H])C([H])([H])[C@]4([H])OC([H])([H])[H])O[C@]3([H])C([H])([H])[H])C([H])([H])[C@@]2([H])OC([H])([H])[H])[C@@]([H])(O[H])[C@]([H])(O[H])[C@@]1([H])O[H]

fn main() {
    let smiles = "[H]C1=C([H])C([H])=C(C([H])([H])N(C(=O)/C([H])=C(\\[H])C(=O)N(N([H])C(=O)[C@@]([H])(N([H])C(=O)[C@@]([H])(N([H])C(=O)OC([H])([H])c2c([H])c([H])c([H])c([H])c2[H])C([H])([H])[H])C([H])([H])[H])C([H])([H])C(=O)N([H])[H])C([H])([H])C2=C([H])C([H])=C([H])O2)O1";
    let mol = ROMol::from_smiles(smiles);
    let leaves = recap_decompose(&mol, None, Some(4), None).get_leaves();
    println!("found {} leaves", leaves.len());
}
