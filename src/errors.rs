use std::fmt::Display;

#[derive(Debug)]
pub struct RDError;

impl Display for RDError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RDError: see stderr for exception info")
    }
}

impl std::error::Error for RDError {}
