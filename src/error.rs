use std::fmt;
pub type FlashResult<T> = Result<T, FlashError>;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FlashError {
    InvalidInput(String),
    OutOfRange(String),
    Degenerate(String),
    NoBracket(String),
    NotConverged(String),
    Json(String),
}

impl fmt::Display for FlashError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FlashError::InvalidInput(message)
            | FlashError::OutOfRange(message)
            | FlashError::Degenerate(message)
            | FlashError::NoBracket(message)
            | FlashError::NotConverged(message)
            | FlashError::Json(message) => f.write_str(message),
        }
    }
}

impl std::error::Error for FlashError {}

impl From<serde_json::Error> for FlashError {
    fn from(error: serde_json::Error) -> Self {
        FlashError::Json(error.to_string())
    }
}
