use std::convert::{From, TryFrom};

use crate::error::ParseError;
use crate::interval::Interval;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct CDS {
    pub range: Interval,
    pub phase: Phase,
}

impl CDS {
    pub fn new(range: Interval, phase: Phase) -> Self {
        Self { range, phase }
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Phase {
    Zero,
    One,
    Two,
}

impl From<Phase> for u8 {
    fn from(phase: Phase) -> u8 {
        match phase {
            Phase::Zero => 0,
            Phase::One => 1,
            Phase::Two => 2,
        }
    }
}

impl From<Phase> for char {
    fn from(phase: Phase) -> Self {
        match phase {
            Phase::Zero => '0',
            Phase::One => '1',
            Phase::Two => '2',
        }
    }
}

const PHASES_STR: &str = "0,1 or 2";

impl TryFrom<&str> for Phase {
    type Error = ParseError;
    fn try_from(phase: &str) -> Result<Self, ParseError> {
        if phase.len() != 1 {
            Err(ParseError::somewhere(PHASES_STR, phase.to_string()))
        } else {
            let c = phase
                .chars()
                .next()
                .expect("string is exactly 1 chars long");
            match c {
                '0' => Ok(Phase::Zero),
                '1' => Ok(Phase::One),
                '2' => Ok(Phase::Two),
                _ => Err(ParseError::somewhere(PHASES_STR, c.to_string())),
            }
        }
    }
}

impl TryFrom<u8> for Phase {
    type Error = ParseError;
    fn try_from(phase: u8) -> Result<Self, ParseError> {
        match phase {
            0 => Ok(Phase::Zero),
            1 => Ok(Phase::One),
            2 => Ok(Phase::Two),
            _ => Err(ParseError::somewhere(PHASES_STR, phase.to_string())),
        }
    }
}

impl TryFrom<char> for Phase {
    type Error = ParseError;
    fn try_from(phase: char) -> Result<Self, ParseError> {
        match phase {
            '0' => Ok(Phase::Zero),
            '1' => Ok(Phase::One),
            '2' => Ok(Phase::Two),
            _ => Err(ParseError::somewhere(PHASES_STR, phase.to_string())),
        }
    }
}

impl std::fmt::Display for Phase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Zero => "0",
                Self::One => "1",
                Self::Two => "2",
            }
        )
    }
}
