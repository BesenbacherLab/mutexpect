use std::fmt;
use std::path::{Path, PathBuf};

use thiserror::Error;

#[derive(Debug, Error)]
#[error("Expected {expected} {location} but observed: {observed}")]
pub struct ParseError {
    expected: &'static str,
    observed: String,
    location: Location,
}

#[derive(Debug)]
pub enum Location {
    Unknown,
    File { path: PathBuf, line: usize },
    Item { type_: &'static str, index: usize },
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Location::Unknown => write!(f, "at unknown location"),
            Location::File { path, line } => {
                write!(f, "in file {} on line {}", path.as_path().display(), line)
            }
            Location::Item { type_, index } => {
                write!(f, "for item of type {} at index {}", type_, index)
            }
        }
    }
}

impl ParseError {
    pub fn somewhere(expected: &'static str, observed: String) -> Self {
        Self {
            expected,
            observed,
            location: Location::Unknown,
        }
    }

    pub fn file(path: PathBuf, line: usize, expected: &'static str, observed: String) -> Self {
        let location = Location::File { path, line };
        Self {
            observed,
            expected,
            location,
        }
    }

    pub fn item(
        type_: &'static str,
        index: usize,
        expected: &'static str,
        observed: String,
    ) -> Self {
        let location = Location::Item { type_, index };
        Self {
            observed,
            expected,
            location,
        }
    }
}

#[derive(Debug, Error)]
pub struct FileError {
    path: Option<PathBuf>,
    #[source]
    source: FileErrorSource,
}

impl FileError {
    pub fn io<P: AsRef<Path>>(path: Option<P>, error: std::io::Error) -> Self {
        let path = match path {
            Some(p) => Some(p.as_ref().to_path_buf()),
            None => None,
        };
        Self {
            path,
            source: error.into(),
        }
    }

    pub fn parse<P: AsRef<Path>>(path: Option<P>, error: ParseError) -> Self {
        let path = match path {
            Some(p) => Some(p.as_ref().to_path_buf()),
            None => None,
        };
        Self {
            path,
            source: error.into(),
        }
    }
}

impl fmt::Display for FileError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.path {
            Some(path) => write!(f, "Failed to work with file {}", path.display()),
            None => write!(f, "Failed to work fith anonymous file"),
        }
    }
}

#[derive(Debug, Error)]
pub enum FileErrorSource {
    #[error("Failed to parse file")]
    Parse {
        #[from]
        source: ParseError,
    },
    #[error("Failed to read/write to file")]
    IO {
        #[from]
        source: std::io::Error,
    },
}

#[derive(Debug, Error)]
pub struct SequenceError {
    sequence_name: String,
    message: &'static str,
    #[source]
    source: Option<Box<dyn std::error::Error>>,
}

impl SequenceError {
    pub fn new(
        sequence_name: String,
        message: &'static str,
        source: Option<Box<dyn std::error::Error>>,
    ) -> Self {
        Self {
            sequence_name,
            message,
            source,
        }
    }
}

impl fmt::Display for SequenceError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Error in Sequence {}: {}",
            self.sequence_name, self.message
        )
    }
}

/// Catch-all error for top-level API
#[derive(Debug, Error)]
pub enum MutexpectError {
    #[error(transparent)]
    ParseError(#[from] ParseError),
    #[error(transparent)]
    FileError(#[from] FileError),
    #[error(transparent)]
    SequenceError(#[from] SequenceError),
}
