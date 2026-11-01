//! The error type of the crate: plain-language messages that name the features involved. Converts into `String` because the
//! femcore APIs use `Result<_, String>`.

use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub struct GeomError(pub String);

pub type GeomResult<T> = Result<T, GeomError>;

impl fmt::Display for GeomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for GeomError {}

impl From<GeomError> for String {
    fn from(e: GeomError) -> String {
        e.0
    }
}

impl From<String> for GeomError {
    fn from(s: String) -> GeomError {
        GeomError(s)
    }
}

impl From<&str> for GeomError {
    fn from(s: &str) -> GeomError {
        GeomError(s.to_string())
    }
}

/// `err!("...", args)` builds an `Err(GeomError)`.
#[macro_export]
macro_rules! err {
    ($($t:tt)*) => { Err($crate::error::GeomError(format!($($t)*))) };
}
