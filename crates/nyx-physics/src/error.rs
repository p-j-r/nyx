//! The error type every `Physics` method returns. Plain data, so it can cross a process boundary later.

use std::fmt;

#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// the problem is malformed or inconsistent; the text is written for the person who made it
    Invalid(String),
    /// the physics does not offer what the problem asks for (a study, a boundary condition, a property)
    Unsupported(String),
    /// a solver failed (singular system, no convergence, ...)
    Solver(String),
    /// a result contained NaN or infinity
    NonFinite(String),
    /// the run was cancelled through `Ctx`
    Cancelled,
    /// the physics crate is a skeleton and its kernels have not been moved in yet (restructure phase 4)
    NotYetMoved(String),
}

impl Error {
    /// A stable short name for the kind, for logs, JSON and MCP error codes.
    pub fn code(&self) -> &'static str {
        match self {
            Error::Invalid(_) => "invalid",
            Error::Unsupported(_) => "unsupported",
            Error::Solver(_) => "solver",
            Error::NonFinite(_) => "non_finite",
            Error::Cancelled => "cancelled",
            Error::NotYetMoved(_) => "not_yet_moved",
        }
    }

    /// Maps the `String` errors of the kernels in `femcore` (the text `"cancelled"` is the cancel signal there).
    pub fn from_kernel(text: String) -> Error {
        if text == "cancelled" {
            Error::Cancelled
        } else if text.contains("NaN or infinity") || text.contains("NaN or infinite") {
            Error::NonFinite(text)
        } else {
            Error::Solver(text)
        }
    }

    pub fn invalid(text: impl Into<String>) -> Error {
        Error::Invalid(text.into())
    }
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Invalid(t)
            | Error::Unsupported(t)
            | Error::Solver(t)
            | Error::NonFinite(t)
            | Error::NotYetMoved(t) => f.write_str(t),
            Error::Cancelled => f.write_str("cancelled"),
        }
    }
}

impl std::error::Error for Error {}
