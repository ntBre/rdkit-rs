fn main() {
    let include = std::env::var("DEP_SHIM_INCLUDE").unwrap();
    println!("cargo:include={include}");
    println!("cargo:rustc-link-arg=-Wl,-rpath,{include}");
}
