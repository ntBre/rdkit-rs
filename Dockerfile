# syntax=docker/dockerfile:1

FROM archlinux AS base
RUN pacman -Syu --noconfirm base-devel boost boost-libs clang freetype2

FROM base AS build-rdkit
RUN <<EOT bash
pacman -Syu --noconfirm git cmake python python-numpy
mkdir /opt
cd /opt
git clone --depth 1 https://github.com/rdkit/rdkit
mkdir rdkit/build
cd rdkit/build
cmake -DRDK_BUILD_INCHI_SUPPORT=ON -DRDK_BUILD_PYTHON_WRAPPERS=OFF -DRDK_BUILD_CPP_TESTS=OFF ..
make -j $(($(nproc) / 2))
EOT

FROM base
COPY --from=build-rdkit /opt/rdkit /opt/rdkit/
RUN pacman -Syu --noconfirm rustup
RUN rustup toolchain install nightly

# Set the entry point to the Rust application
CMD ["/app/target/release/rdkit-rs"]
