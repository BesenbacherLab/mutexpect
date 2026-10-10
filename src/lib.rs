mod amino_acid;
mod cds;
pub mod error;
pub mod interval;
mod mutation;
mod point_mutation_classifier;
mod seq_window_slider;
mod sequence_annotation;

use std::collections::hash_map::HashMap;
use std::convert::{From, TryFrom};
use std::fmt;

use serde_repr::{Deserialize_repr, Serialize_repr};

use pattern_partition_prediction::{PaPaPred, PaPaPredIndel};
use twobit::TwoBitFile;

pub use crate::cds::{Phase, CDS};
use crate::error::{MutexpectError, ParseError};
pub use crate::interval::Interval;
use crate::mutation::{PointMutation, Indel};
pub use crate::point_mutation_classifier::PointMutationClassifier;
use crate::seq_window_slider::SeqWindowSlider;
pub use crate::sequence_annotation::{read_sequence_annotations_from_file, SeqAnnotation, Strand};

/// What a possible mutation changes at a position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Alteration {
    /// A point mutation to this base.
    Point(Base),
    /// An in-frame insertion or deletion (anchored at this position).
    InFrameIndel,
    /// A frameshift insertion or deletion (anchored at this position).
    FrameshiftIndel,
}

/// A possible mutation and its class, see [`for_each_possible_mutation`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PossibleMutation {
    pub alteration: Alteration,
    pub mutation_type: MutationType,
}

/// One position of an annotated sequence and the mutations possible there.
#[derive(Debug)]
pub struct Site<'a> {
    /// 0-based genomic position.
    pub genomic_position: usize,
    pub ref_base: Base,
    /// The `2 * seq_flanking_length + 1` bases centred on the position.
    pub seq_window: &'a str,
    /// Point mutations (to each other base, in A, C, G, T order) that pass the
    /// filters, followed by indels (only inside coding sequences).
    pub mutations: &'a [PossibleMutation],
}

/// Call `f` for every position of `seq_annotation` with the mutations that are
/// possible there and their classes.
///
/// `seq` is the genomic sequence of the annotated region plus
/// `seq_flanking_length` bases on each side. `f` is called for every position,
/// also when no mutation passes the filters (then `mutations` is empty).
/// Intronic and unknown point mutations are only included if requested. With
/// `filter_plof`, frameshift indels in the last 50 bases of the coding sequence
/// are left out (and coding point mutations classified accordingly).
pub fn for_each_possible_mutation<E>(
    seq: &str,
    seq_annotation: &SeqAnnotation,
    seq_flanking_length: usize,
    include_intronic: bool,
    include_unknown: bool,
    filter_plof: bool,
    mut f: impl FnMut(&Site) -> Result<(), E>,
) -> Result<(), E> {
    let window_length = 1 + 2 * seq_flanking_length;
    let sequence_genomic_start_position = seq_annotation.range.start;
    let classifier = PointMutationClassifier::new(&seq_annotation, seq_flanking_length);
    let introns = seq_annotation.get_introns();
    let mut introns_iter = introns.iter();
    let mut cds_iter = seq_annotation.coding_sequences.iter();

    let mut current_intron: Option<Interval> = introns_iter.next().copied();
    let mut current_cds: Option<CDS> = cds_iter.next().copied();
    let mut mutations = Vec::with_capacity(5);
    for (i, seq_window) in SeqWindowSlider::new(seq, window_length)
        .into_iter()
        .enumerate()
    {
        let genomic_position = sequence_genomic_start_position + i;
        let seq_window_vec: Vec<char> = seq_window.chars().collect();
        let ref_base: Base = seq_window_vec[seq_flanking_length].into();

        let overlapping_intron = {
            loop {
                if let Some(intron) = current_intron {
                    if genomic_position < intron.start {
                        break None; // still has some way to go to the next intron
                    } else if genomic_position < intron.stop {
                        break Some(intron); // yay, we are inside the current intron!
                    } else {
                        current_intron = introns_iter.next().copied();
                    }
                } else {
                    break None; // iterator is exhausted
                }
            }
        };

        let overlapping_cds = {
            loop {
                if let Some(cds) = &current_cds {
                    if genomic_position < cds.range.start {
                        break None;
                    } else if genomic_position < cds.range.stop {
                        break Some(cds); // yay, we are inside the current cds!
                    } else {
                        current_cds = cds_iter.next().copied();
                    }
                } else {
                    break None; // iterator is exhausted
                }
            }
        };

        let mut mutation_type = classifier.classify_by_position(
            genomic_position,
            &seq_window_vec,
            &overlapping_intron, // may be none
        );

        mutations.clear();
        for other_nuc in &[Base::A, Base::C, Base::G, Base::T] {
            if *other_nuc == ref_base {
                continue;
            }
            if let Some(cds) = overlapping_cds {
                // if we are in a coding region we need to do additional classifications
                mutation_type = classifier.classify_coding_mutation(
                    genomic_position,
                    &seq_window_vec,
                    other_nuc.name(),
                    &cds,
                    filter_plof,
                )
            } // else just recycle the mutation class
            if mutation_type == MutationType::Intronic && !include_intronic {
                continue;
            }
            if mutation_type == MutationType::Unknown && !include_unknown {
                continue;
            }
            mutations.push(PossibleMutation {
                alteration: Alteration::Point(*other_nuc),
                mutation_type,
            });
        }

        // indels: we don't care about intronic indels because most of them don't matter
        if let Some(cds) = overlapping_cds {
            if cds.range.start < genomic_position {
                // needs to be greater, because of anchor base
                mutations.push(PossibleMutation {
                    alteration: Alteration::InFrameIndel,
                    mutation_type: MutationType::InFrameIndel,
                });
                let in_last_50 = match seq_annotation.strand {
                    Strand::Plus => genomic_position >= seq_annotation.get_end_minus50(),
                    Strand::Minus => genomic_position <= seq_annotation.get_end_minus50(),
                };
                if !(filter_plof && in_last_50) {
                    mutations.push(PossibleMutation {
                        alteration: Alteration::FrameshiftIndel,
                        mutation_type: MutationType::FrameshiftIndel,
                    });
                }
            }
        }

        f(&Site {
            genomic_position,
            ref_base,
            seq_window,
            mutations: &mutations,
        })?;
    }
    Ok(())
}

pub fn possible_mutations(
    seq: &str,
    seq_annotation: &SeqAnnotation,
    point_mutation_probabilities: &PaPaPred,
    indel_probabilities: &Option<PaPaPredIndel>,
    drop_nan: bool,
    include_intronic : bool,
    include_unknown : bool,
    filter_plof : bool,
) -> Result<Vec<MutationEvent>, MutexpectError> {
    let mut result = Vec::new();
    let seq_flanking_length = point_mutation_probabilities.kmer_size() / 2; // this should round down
    for_each_possible_mutation(
        seq,
        seq_annotation,
        seq_flanking_length,
        include_intronic,
        include_unknown,
        filter_plof,
        |site| -> Result<(), MutexpectError> {
            let mutation_probability = point_mutation_probabilities
                .rates(site.seq_window)
                .map_err(|e| {
                    error::SequenceError::new(
                        seq_annotation.name.clone(),
                        "Bad sequence (for point mutation)",
                        Some(Box::new(e)),
                    )
                })?;
            let mut indel_rates = None;
            for m in site.mutations {
                let probability = match m.alteration {
                    Alteration::Point(other_nuc) => {
                        let probability = mutation_probability.by_name(other_nuc.name());
                        if drop_nan && probability.is_nan() {
                            // avoid trouble but make it impossible to reconstruct the positions of the point mutations
                            // but that's not something I'm concerned with at the moment anyway.
                            continue;
                        }
                        probability
                    }
                    Alteration::InFrameIndel | Alteration::FrameshiftIndel => {
                        let p = match indel_probabilities {
                            Some(p) => p,
                            None => continue,
                        };
                        if indel_rates.is_none() {
                            let w = site.seq_window;
                            indel_rates = Some(p.rates(&w[0..w.len() - 1]) // 01[2]34 -> 01[]23
                                .map_err(|e| {
                                    error::SequenceError::new(
                                        seq_annotation.name.clone(),
                                        "Bad sequence (for indel split)",
                                        Some(Box::new(e)),
                                    )
                                })?);
                        }
                        let rates = indel_rates.as_ref().expect("set above");
                        if m.alteration == Alteration::InFrameIndel {
                            rates.inframe
                        } else {
                            rates.outframe
                        }
                    }
                };
                result.push(MutationEvent {
                    mutation_type: m.mutation_type,
                    probability,
                });
            }
            Ok(())
        },
    )?;
    Ok(result)
}

pub fn expected_number_of_mutations(
    possible_mutations: &[MutationEvent],
) -> HashMap<MutationType, f64> {
    let mut result = HashMap::new();
    for event in possible_mutations {
        *result.entry(event.mutation_type).or_insert(0.0) += event.probability as f64
    }
    result
}

pub fn observed_number_of_mutations(
    mutations: &[PointMutation], // point mutations within the seq_annotation (the gene/transcript)
    indels: &[Indel],
    seq_annotation: &SeqAnnotation, //only one annotation for the current gene/transcript
    twobit_ref_seq: &TwoBitFile,
    filter_plof : bool,
) -> Result<HashMap<MutationType, usize>, MutexpectError> {
    let mut result = HashMap::new();
    let dna = twobit_ref_seq; // alias

    // we will provide a 5 nucleotide window around the mutation
    let classifier = PointMutationClassifier::new(seq_annotation, 5);

    for mutation in mutations {
        let start = mutation.position - 2; // 2 nucs flanking to the left
        let stop = mutation.position + 3; // 2 nucs flanking to the right
        let nucleotides: Vec<char> = dna
            .sequence(&mutation.chromosome, start, stop)
            .unwrap()
            .chars()
            .collect();
        let mut mut_type = classifier.classify_by_position(
            mutation.position,
            &nucleotides,
            &seq_annotation.find_intron(mutation.position),
        );
        if mut_type == MutationType::Unknown {
            if let Some(cds) = seq_annotation.find_cds(mutation.position) {
                mut_type = classifier.classify_coding_mutation(
                    mutation.position,
                    &nucleotides,
                    mutation.alternative,
                    &cds,
                    filter_plof,
                );
            } // else: not coding
        } // else: just recycle Unknown 3 times
        result.entry(mut_type).and_modify(|c| *c += 1).or_insert(1);
    }

    for indel in indels {
        let mut mut_type = MutationType::Unknown;
        for intron in seq_annotation.get_introns() {
            if intron.contains(indel.position) {
                mut_type = MutationType::Intronic;
                break;
            }
        }

        if mut_type == MutationType::Unknown { // it's not intronic. It might be a frameshift
            for cds in &seq_annotation.coding_sequences {
                if cds.range.contains(indel.position) {
                    if indel.is_inframe() {
                        mut_type = MutationType::InFrameIndel
                    } else {
                        mut_type = MutationType::FrameshiftIndel;
                        if filter_plof {
                            match seq_annotation.strand {
                                Strand::Plus => {
                                    if indel.position >= seq_annotation.get_end_minus50() {
                                        mut_type = MutationType::Unknown
                                    }  
                                },      
                                Strand::Minus => {
                                    if indel.position <= seq_annotation.get_end_minus50() {
                                        mut_type = MutationType::Unknown
                                    } 
                                },
                            }
                        }
                    }
                    break
                }
            }
        }
        result.entry(mut_type).and_modify(|c| *c += 1).or_insert(1);
    }
    Ok(result)
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum MutationType {
    Unknown = 0,
    Synonymous = 1,
    Missense = 2,
    Nonsense = 3,
    StopLoss = 4,
    StartCodon = 5,
    SpliceSite = 6,
    Intronic = 7,
    InFrameIndel = 8,
    FrameshiftIndel = 9,
}

impl MutationType {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Synonymous => "synonymous",
            Self::Missense => "missense",
            Self::Nonsense => "nonsense",
            Self::StopLoss => "stop_loss",
            Self::StartCodon => "start_codon",
            Self::SpliceSite => "splice_site",
            Self::Intronic => "intronic",
            Self::InFrameIndel => "in_frame_indel",
            Self::FrameshiftIndel => "frameshift_indel",
        }
    }

    pub fn iter() -> MutationTypeIter {
        MutationTypeIter { index: 0 }
    }
}


impl fmt::Display for MutationType {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // I get a stack overflow if I simply do self.to_string(). So I'll implement the string
        // representation directly instead.
        // I need to figure out why that happens though... Seems like it's related to different
        // formatting functions calling each other recursively.
        let string = match self {
            Self::Unknown => "Unknown",
            Self::Synonymous => "Synonymous",
            Self::Missense => "Missense",
            Self::Nonsense => "Nonsense",
            Self::StopLoss => "StopLoss",
            Self::StartCodon => "StartCodon",
            Self::SpliceSite => "SpliceSite",
            Self::Intronic => "Intronic",
            Self::InFrameIndel => "InFrameIndel",
            Self::FrameshiftIndel => "FrameshiftIndel",
        };
        write!(f, "{}", string)
    }
}

impl TryFrom<&str> for MutationType {
    type Error = ParseError;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        Ok(match s.to_lowercase().as_str() {
            "unknown" => Self::Unknown,
            "synonymous" => Self::Synonymous,
            "missense" => Self::Missense,
            "nonsense" => Self::Nonsense,
            "stoploss" | "stop_loss" => Self::StopLoss,
            "startcodon" | "start_codon" => Self::StartCodon,
            "splicesite" | "splice_site" => Self::SpliceSite,
            "intronic" => Self::Intronic,
            "in_frame" | "in_frame_indel" | "inframe" => Self::InFrameIndel,
            "out_frame" | "out_frame_outdel" | "outframe" => Self::FrameshiftIndel,
            _ => {
                return Err(ParseError::somewhere(
                    "name of mutation type",
                    s.to_string(),
                ))
            }
        })
    }
}

impl From<u8> for MutationType {
    fn from(n: u8) -> Self {
        match n {
            0 => Self::Unknown,
            1 => Self::Synonymous,
            2 => Self::Missense,
            3 => Self::Nonsense,
            4 => Self::StopLoss,
            5 => Self::StartCodon,
            6 => Self::SpliceSite,
            7 => Self::Intronic,
            8 => Self::InFrameIndel,
            9 => Self::FrameshiftIndel,
            _ => Self::Unknown,
        }
    }
}

pub struct MutationTypeIter {
    index: u8,
}

impl std::iter::Iterator for MutationTypeIter {
    type Item = MutationType;

    fn next(&mut self) -> Option<Self::Item> {
        let mut_type: MutationType = self.index.into();
        self.index += 1;
        if mut_type == MutationType::Unknown && self.index != 1 { // 1, because we already incremented
            None
        } else {
            Some(mut_type)
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MutationEvent {
    pub mutation_type: MutationType,
    pub probability: f32,
}

impl MutationEvent {
    pub fn new(mutation_type: MutationType, probability: f32) -> Self {
        Self {
            mutation_type,
            probability,
        }
    }
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Base {
    pub fn name(&self) -> char {
        match self {
            Base::A => 'A',
            Base::C => 'C',
            Base::G => 'G',
            Base::T => 'T',
            Base::N => 'N',
        }
    }
}

impl From<u8> for Base {
    fn from(byte: u8) -> Self {
        match byte {
            65 | 97 => Self::A,
            67 | 99 => Self::C,
            71 | 103 => Self::G,
            84 | 116 => Self::T,
            78 | 110 => Self::N,
            _ => panic!("Bad nucleotide: {}", byte),
        }
    }
}

impl From<char> for Base {
    fn from(c: char) -> Self {
        match c {
            'A' | 'a' => Self::A,
            'C' | 'c' => Self::C,
            'G' | 'g' => Self::G,
            'T' | 't' => Self::T,
            'N' | 'n' => Self::N,
            _ => panic!("Bad nucleotide: {}", c),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mutation_type_iterator() {
        let mutation_types = vec![
            MutationType::Unknown,
            MutationType::Synonymous,
            MutationType::Missense,
            MutationType::Nonsense,
            MutationType::StopLoss,
            MutationType::StartCodon,
            MutationType::SpliceSite,
            MutationType::Intronic,
            MutationType::InFrameIndel,
            MutationType::FrameshiftIndel,
        ];
        let mut iter = MutationType::iter();
        for mut_type in mutation_types {
            assert_eq!(mut_type, iter.next().unwrap());
        }
        assert!(iter.next().is_none());
    }
}
