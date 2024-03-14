fn main() {
    let include = std::env::var("DEP_SHIM_INCLUDE").unwrap();
    println!("cargo:include={include}")
}
