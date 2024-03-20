[![test](https://github.com/ntBre/rdkit-rs/actions/workflows/test.yml/badge.svg)](https://github.com/ntBre/rdkit-rs/actions/workflows/test.yml)

# rdkit-rs
[RDKit][rdkit] interface in safe Rust

## Usage

### Obtaining RDKit

I have not been able to get this to work with an RDKit installed with conda, so
you must currently build RDKit yourself. To build from source, run this
somewhere on your computer:

``` shell
git clone --depth 1 https://github.com/rdkit/rdkit
cd rdkit
mkdir build
cd build
cmake .. -DRDK_BUILD_INCHI_SUPPORT=ON
make # optionally with -j and however many CPUs you want to give it
```

`rdkit-sys`, and thus `rdkit-rs`, includes InChi key support, so you must
compile RDKit with InChi support enabled.

### Running dependent crates

The `rdkit-sys` build script should take care of building the shared library and
including the absolute path in the binary itself via `rpath` linker arguments.
However, it does require the environment variable `RDROOT` to be set to locate
the RDKit libraries and header files from the steps above. Probably the easiest
way to do this is to add a Cargo config file in `.cargo/config.toml` in your
project:

``` toml
[env]
RDROOT = "/path/to/your/rdkit"
```

But you can also prefix your Cargo commands with `RDROOT=/path/to/rdkit` or
however else you like to set environment variables.

### Complete Example

Initialize a new Rust project and enter the directory:

``` shell
cargo new rdkit-test
cd rdkit-test
```

Add `rdkit-rs` as a dependency:

``` shell
cargo add --git=https://github.com/ntBre/rdkit-rs
```

Place the path to your cloned RDKit code in `.cargo/config.toml`:

``` toml
[env]
RDROOT = "/home/brent/omsf/clone/rdkit"
```

Here I have used my actual path, but be sure to replace it with your own.

Update `src/main.rs`:

``` shell
cat > src/main.rs <<INP
use rdkit_rs::ROMol;

fn main() {
    let mol = ROMol::from_smiles("CCO");
    let inchi = mol.to_inchi_key();
    println!("Hello inchi: {inchi}!")
}
INP
```

Run with cargo:

``` shell
cargo run
```

### Using Docker

You can use the included `Dockerfile` to build a base Arch Linux image with
RDKit and Rust:

``` shell
docker build -t rdkit-rs - < Dockerfile
```

Then you can run the image

``` shell
docker run -it rdkit-rs bash
```

And follow the rest of the steps in the example above. The Dockerfile copies
RDKit to `/opt/rdkit`, so you can use `RDROOT=/opt/rdkit` in your
`.cargo/config.toml` or at the command line.

[rdkit]: https://github.com/rdkit/rdkit
