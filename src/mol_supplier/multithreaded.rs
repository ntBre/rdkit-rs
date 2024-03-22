use std::ffi::CString;

use crate::ROMol;

/// WARNING: use this at your own risk. I don't think I did anything wrong in
/// the C++ wrapper, and it's using the default of 2 writer threads, but I've
/// gotten segfaults reading the ChEMBL 33 SDF. RDKit should have been built
/// with the default RDK_BUILD_THREADSAFE_SSS=ON, so defining that macro in the
/// wrapper should be safe
pub struct MultithreadedSDMolSupplier(
    *mut rdkit_sys::RDKit_MultithreadedSDMolSupplier,
);

unsafe impl Send for MultithreadedSDMolSupplier {}

impl MultithreadedSDMolSupplier {
    /// construct an [MultithreadedSDMolSupplier] from a filepath that can be
    /// converted to a CString. panics if this conversion fails
    pub fn new(path: impl Into<Vec<u8>>) -> Self {
        let cpath = CString::new(path).expect("failed to create CString");
        unsafe {
            let inner =
                rdkit_sys::RDKit_MultithreadedSDMolSupplier_new(cpath.as_ptr());
            Self(inner)
        }
    }

    /// reports whether or not `self` is at the end of the underlying file
    pub fn at_end(&self) -> bool {
        unsafe { rdkit_sys::RDKit_MultithreadedSDMolSupplier_at_end(self.0) }
    }
}

impl Drop for MultithreadedSDMolSupplier {
    fn drop(&mut self) {
        unsafe {
            rdkit_sys::RDKit_MultithreadedSDMolSupplier_delete(self.0);
        }
    }
}

impl Iterator for MultithreadedSDMolSupplier {
    type Item = ROMol;

    fn next(&mut self) -> Option<Self::Item> {
        if self.at_end() {
            return None;
        }
        Some(unsafe {
            ROMol(rdkit_sys::RDKit_MultithreadedSDMolSupplier_next(self.0))
        })
    }
}
