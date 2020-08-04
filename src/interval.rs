use std::cmp::{max, min};
use std::fmt;

use crate::error::ParseError;

/// A range of values
///
/// By convention, start and stop are zero-based.
/// The start position is always inclusive and the stop position is always exclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Interval {
    pub start: usize, // inclusive
    pub stop: usize,  // exclusive
}

impl Interval {
    pub fn new(start: usize, stop: usize) -> Result<Interval, ParseError> {
        if stop < start {
            Err(ParseError::somewhere(
                "low < high",
                format!("{}>={}", start, stop),
            ))
        } else {
            Ok(Interval { start, stop })
        }
    }

    pub fn parse(range: &str) -> Result<Interval, ParseError> {
        let split = range.split('-');
        let parts: Vec<&str> = split.collect();
        if parts.len() != 2 {
            return Err(ParseError::somewhere("start-stop", range.to_string()));
        }
        // we now know for sure that parts is of length 2

        if let Ok(start) = parts[0].parse::<usize>() {
            if let Ok(stop) = parts[1].parse::<usize>() {
                if stop <= start {
                    Err(ParseError::somewhere("start<stop", range.to_string()))
                } else {
                    Interval::new(start, stop)
                }
            } else {
                Err(ParseError::somewhere("usize", parts[1].to_string()))
            }
        } else {
            Err(ParseError::somewhere("usize", parts[0].to_string()))
        }
    }

    pub fn parse_many(text: &str, delimiter: char) -> Result<Vec<Interval>, ParseError> {
        let mut result = Vec::new();
        if text.trim() != "" {
            // if the input is an empty string, return an empty vector
            for token in text.split(delimiter) {
                result.push(Interval::parse(token)?)
            }
        }
        Ok(result)
    }

    pub fn contains(&self, pos: usize) -> bool {
        self.start <= pos && pos < self.stop
    }

    pub fn intersection(&self, other: &Interval) -> Option<Interval> {
        let start = max(self.start, other.start);
        let stop = min(self.stop, other.stop);
        if start == stop {
            None // interval of length 0 is as good as no interval at all
        } else {
            Interval::new(start, stop).ok() // this will be None, if there is no intersection
        }
    }

    pub fn len(&self) -> usize {
        self.stop - self.start
    }

    pub fn is_empty(&self) -> bool {
        self.stop == self.start
    }
}

impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.start, self.stop)
    }
}
