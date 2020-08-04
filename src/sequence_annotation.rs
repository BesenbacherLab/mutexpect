use std::convert::{TryFrom, TryInto};
use std::path::Path;

use tabfile::Tabfile;

use crate::cds::CDS;
use crate::error::{FileError, ParseError};
use crate::interval::Interval;

#[derive(Debug, PartialEq)]
pub struct SeqAnnotation {
    pub name: String,               // col 0
    pub chr: String,                // col 1
    pub range: Interval,            // col 2
    pub strand: Strand,             // col 3
    pub exons: Vec<Interval>,       // col 4
    pub coding_sequences: Vec<CDS>, // col 5+6
}

impl SeqAnnotation {
    pub fn new(
        name: String,
        chr: String,
        range: Interval,
        strand: Strand,
        exons: Vec<Interval>,
        coding_sequences: Vec<CDS>,
    ) -> Self {
        Self {
            name,
            chr,
            range,
            strand,
            exons,
            coding_sequences,
        }
    }

    /// Find introns based on the exons
    /// Assumption: Every intron is flanked by exons
    ///
    /// Cases:  a    b    c    d    e
    /// exons   |  --|--  |  --|--  |
    /// DNA  --------------------------
    ///
    /// a) Not an intron, because it is before the first exon
    /// b,d) Not an intron, because it is inside an exon
    /// c) An intron, because it is between two exons
    /// e) 3' UTR region. Not an intron.
    ///
    pub fn find_intron(&self, genomic_position: usize) -> Option<Interval> {
        let mut left_exon_end = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            if genomic_position < exon.start {
                if i > 0 {
                    // case c
                    return Some(Interval::new(left_exon_end, exon.start).expect("start<stop"));
                } else {
                    // case a
                    return None;
                }
            } else if exon.contains(genomic_position) {
                return None; // cases b and d
            } // the position is right of the current exon
            left_exon_end = exon.stop;
        }
        None // case e
    }

    /// Determine all intron locations based on the exons
    ///
    /// Note that exons may overlap. But in any case, this function only returns
    /// regions that are not covered by any exons (assuming that the exons are
    /// sorted by start and then by stop position).
    pub fn get_introns(&self) -> Vec<Interval> {
        let mut result = Vec::new();
        let mut left_exon_end = 0;
        for (i, exon) in self.exons.iter().enumerate() {
            if i > 0 && left_exon_end < exon.start {
                // exons do not overlap. We can define a pure intron region
                result.push(
                    Interval::new(left_exon_end, exon.start).unwrap(), // pretty sure the exons are ordered
                )
            } // else: overlapping exons. We need to treat this as a super-exon region
            left_exon_end = exon.stop
        }
        result
    }

    pub fn find_exon(&self, genomic_position: usize) -> Option<Interval> {
        for exon in &self.exons {
            if exon.contains(genomic_position) {
                return Some(*exon);
            }
        }
        None
    }

    pub fn find_cds(&self, genomic_position: usize) -> Option<CDS> {
        for cds in &self.coding_sequences {
            if cds.range.contains(genomic_position) {
                return Some(*cds);
            }
        }
        None
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Plus => "+",
                Self::Minus => "-",
            }
        )
    }
}

impl TryFrom<char> for Strand {
    type Error = ParseError;
    fn try_from(c: char) -> Result<Self, Self::Error> {
        match c {
            '+' => Ok(Strand::Plus),
            '-' => Ok(Strand::Minus),
            _ => Err(ParseError::somewhere("+ or -", c.to_string())),
        }
    }
}

pub fn read_sequence_annotations_from_file<P: AsRef<Path>>(
    path: P,
    filter_for_id: Option<&str>,
) -> Result<Vec<SeqAnnotation>, FileError> {
    const NAME_IDX: usize = 0;
    const CHR_IDX: usize = 1;
    const STRAND_IDX: usize = 2;
    const GENE_RANGE_IDX: usize = 3;
    const EXONS_IDX: usize = 4;
    const CDS_IDX: usize = 5;
    const CDS_PHASES_IDX: usize = 6;

    let mut result = Vec::new();
    let tabfile = match Tabfile::open(&path) {
        Ok(tf) => tf.comment_character('#'),
        Err(e) => return Err(FileError::io(Some(&path), e)),
    };
    for record_result in tabfile {
        let record = match record_result {
            Ok(record) => record,
            Err(e) => return Err(FileError::io(Some(&path), e)),
        };
        let tokens = record.fields();
        if tokens.len() < 7 {
            let err = ParseError::file(
                path.as_ref().to_path_buf(),
                record.line_number(),
                "7 columns",
                record.line().to_string(),
            );
            return Err(FileError::parse(Some(&path), err));
        }

        if let Some(id) = filter_for_id {
            if tokens[NAME_IDX] != id {
                continue;
            }
        }
        let exons = Interval::parse_many(tokens[EXONS_IDX], ';')
            .map_err(|e| FileError::parse(Some(&path), e))?;
        let cds_intervals = Interval::parse_many(tokens[CDS_IDX], ';')
            .map_err(|e| FileError::parse(Some(&path), e))?;
        let cds_phases: Vec<&str> = tokens[CDS_PHASES_IDX]
            .split(';')
            .filter(|p| p != &"")
            .collect();
        if cds_intervals.len() != cds_phases.len() {
            return Err(FileError::parse(
                Some(&path),
                ParseError::file(
                    path.as_ref().to_path_buf(),
                    record.line_number(),
                    "number of CDS regions and CDS phases to be the same",
                    format!(
                        "#CDS_regions={}, #CDS_phases={}",
                        cds_intervals.len(),
                        cds_phases.len()
                    ),
                ),
            ));
        }
        let coding_sequences = {
            let mut cdss = Vec::new();
            for (interval, phase) in cds_intervals.iter().zip(cds_phases) {
                let phase = match phase.try_into() {
                    Ok(phase) => phase,
                    Err(e) => return Err(FileError::parse(Some(&path), e)),
                };
                cdss.push(CDS::new(*interval, phase));
            }
            cdss
        };
        let strand: char = tokens[STRAND_IDX].chars().next().unwrap();
        result.push(SeqAnnotation {
            name: tokens[NAME_IDX].to_string(),
            chr: tokens[CHR_IDX].to_string(),
            range: Interval::parse(tokens[GENE_RANGE_IDX])
                .map_err(|e| FileError::parse(Some(&path), e))?,
            strand: strand
                .try_into()
                .map_err(|e| FileError::parse(Some(&path), e))?,
            exons,
            coding_sequences,
        });
    }
    Ok(result)
}
