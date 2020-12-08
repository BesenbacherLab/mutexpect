#![allow(clippy::identity_op)]

use std::fmt;

use crate::amino_acid::{translate, AminoAcid};
use crate::cds::CDS;
use crate::interval::Interval;
use crate::sequence_annotation::{SeqAnnotation, Strand};
use crate::MutationType;

/// Ensure that modulo also works correctly for negative numbers
fn mod3(i: isize) -> usize {
    i.rem_euclid(3) as usize
}

fn complement(nuc: char) -> char {
    match nuc {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        _ => panic!("Unsupported nucleotide {}", nuc),
    }
}

pub struct PointMutationClassifier<'a> {
    gene_info: &'a SeqAnnotation,
    middle_index: usize,
}

impl<'a> PointMutationClassifier<'a> {
    pub fn new(gene_info: &'a SeqAnnotation, middle_index: usize) -> Self {
        Self {
            gene_info,
            middle_index,
        }
    }

    /// Classify the effect of a point mutation
    ///
    /// Parameters:
    /// `genomic_position`: The orginal genomic position on a chromosome where the mutation takes place
    /// `seq`: The reference sequence around the mutation site, with the mutated locus in the middle
    /// `other_nucleotide`: The nucleotide that will replace the old nucleotide
    /// `cds_index_hint`: (optional) The internal array-index of the CDS that this mutation falls into. If you provide this, and it is correct, no search for a matching CDS needs to be done.
    ///
    /// Returns: The corresponding `MutationType`
    pub fn classify_by_position(
        &self,
        genomic_position: usize,
        seq: &[char],
        intron: &Option<Interval>,
    ) -> MutationType {
        if intron.is_some() {
            if self.disrupts_cannonical_splice_site(genomic_position, seq, &intron.expect("is some")) {
                MutationType::SpliceSite
            } else {
                MutationType::Intronic
            }
        } else if self.disrupts_start_codon(genomic_position, seq) {
            MutationType::StartCodon
        } else {
            MutationType::Unknown
        }
    }

    /// Returns `true` if the point mutation disrupts a cannonical splice site.
    /// A cannonical splice site is a left "GT" or a right "AG" at the ends of
    /// an intron.
    ///
    /// There are basically 4 locations of interest for this test:
    ///
    ///    12    34
    /// ...GT....AG...
    ///    ^      ^
    ///    |intron|
    ///
    /// If any of these four bases are mutated, return True. Also take into
    /// account the reverse complement:
    ///
    ///    |intron|
    ///    v      v
    /// ...GT....AG...
    /// ...CA....TC...
    ///    43    21 <- reverse complement tests case indexes
    ///
    fn disrupts_cannonical_splice_site(
        &self,
        genomic_position: usize,
        seq: &[char],
        intron: &Interval,
    ) -> bool {
        let middle_index = self.middle_index;
        let start = intron.start;
        let stop = intron.stop;

        match self.gene_info.strand {
            Strand::Plus => {
                // 5' to 3' direction
                if genomic_position == start {
                    // check for [G]T
                    seq[middle_index] == 'G' && seq[middle_index + 1] == 'T'
                } else if genomic_position == start + 1 {
                    // check for G[T]
                    seq[middle_index - 1] == 'G' && seq[middle_index] == 'T'
                } else if genomic_position == stop - 2 {
                    // check for [A]G
                    seq[middle_index] == 'A' && seq[middle_index + 1] == 'G'
                } else if genomic_position == stop - 1 {
                    // check for A[G]
                    seq[middle_index - 1] == 'A' && seq[middle_index] == 'G'
                } else {
                    false // not a canonical splice site
                }
            }
            Strand::Minus => {
                // Because the reference sequence is on the plus strand, we need to check for the reverse complement
                if genomic_position == start {
                    // check for [C]T
                    seq[middle_index] == 'C' && seq[middle_index + 1] == 'T'
                } else if genomic_position == start + 1 {
                    // check for C[T]
                    seq[middle_index - 1] == 'C' && seq[middle_index - 0] == 'T'
                } else if genomic_position == stop - 2 {
                    // check for [A]C
                    seq[middle_index] == 'A' && seq[middle_index + 1] == 'C'
                } else if genomic_position == stop - 1 {
                    // check for A[C]
                    seq[middle_index - 1] == 'A' && seq[middle_index - 0] == 'C'
                } else {
                    false // not a canonical splice site
                }
            }
        }
    }

    /// Check if a mutation at this position would disrupt the first codon
    ///
    /// Parameters:
    ///   genomic_position: The position on the chromosome
    fn disrupts_start_codon(&self, genomic_position: usize, _sequence: &[char]) -> bool {
        let gene = self.gene_info;
        match gene.strand {
            Strand::Plus => {
                if let Some(cds) = gene.coding_sequences.first() {
                    // the first CDS usually starts with the start codon
                    if let Some(distance_to_start_codon) =
                        genomic_position.checked_sub(cds.range.start)
                    {
                        /* //too much output for now
                        if distance_to_start_codon == 0 { // then we do a sanity check
                            let flanking: usize = _sequence.len() / 2;
                            if _sequence[ flanking .. flanking + 3 ] !=  [ 'A', 'T', 'G' ] {
                                eprintln!( "[WARNING] Lowest CDS on the + strand does not have a proper start codon. This can be correct for some genes" );
                            }
                        }
                        */
                        return distance_to_start_codon < 3; // usize is always >= 0
                    }
                }
            }
            Strand::Minus => {
                if let Some(cds) = gene.coding_sequences.last() {
                    // the last CDS must end with a reverse start codon
                    if let Some(distance_to_start_codon) =
                        cds.range.stop.checked_sub(genomic_position)
                    {
                        /* //too much output for now
                        if distance_to_start_codon == 3 {
                            let flanking: usize = _sequence.len() / 2;
                            if _sequence[ flanking .. flanking + 3 ] != [ 'C', 'A', 'T' ] {
                                eprintln!( "[WARNING] Highest CDS on the - strand does not have a proper start codon. This can be correct for some genes" );
                            }
                        }
                        */
                        return 0 < distance_to_start_codon && distance_to_start_codon < 4;
                        // keep in mind that the .stop is exclusive
                    }
                }
            }
        }
        false
    }

    pub fn classify_coding_mutation(
        &self,
        genomic_position: usize,
        seq: &[char],
        other_nucleotide: char,
        cds: &CDS,
    ) -> MutationType {
        let phase = cds.phase as usize;
        let middle_index = self.middle_index;

        let (old_codon, relative_frame) = match self.gene_info.strand {
            Strand::Plus => {
                // easy mode
                // The relative frame is the coding frame of the `genomic_position`
                let relative_frame = mod3(
                    // this is ugly.
                    genomic_position as isize - cds.range.start as isize - phase as isize,
                ); // I need to establish a baseline for the frames. That baseline is (cds.range.start + phase). Then I need to calculate the distance from the baseline to the mutation

                debug_assert!(relative_frame < 3, "relative_frame = {}", relative_frame);
                let pos_of_codon_start = middle_index - relative_frame;
                let nucs: String = seq[pos_of_codon_start..pos_of_codon_start + 3]
                    .iter()
                    .collect();
                debug_assert_eq!(nucs.len(), 3);
                (nucs, relative_frame)
            }
            Strand::Minus => {
                // upside down
                //I want to have the following properties for the frame:
                //relative_frame=0 => first base of codon
                //relative_frame=1 => second base of codon
                //relative_frame=2 => third base of codon
                //
                //This means that the frames are as follows along the sequence:
                //
                // plus strand:  5'     non-coding     3'
                // minus strand: 3'   210 210 210 210  5'
                //
                // And we also process the CDS from the 5' to the 3' on the minus strand
                //
                // Example for a frame calculation:
                //
                // [210][210][2$X]   reverse strand frames ($=CDS end position)
                //  123  456  789    genomic positons
                //            ^
                //            |--- last inclusive position. So we have a phase of 1.
                //
                // Let's say we are at position 4. Then we need to know the
                // distance to the end of the cds which is (8-1)-4=3. Then we
                // have to take into account the phase which in this case is 1
                // at $, so we have 3-1 = 2. And 2 mod 3 is 2.
                //
                // Another example:
                //
                // [210][21$]
                //  123  456
                //        ^
                //        |--- last inclusive position. We have a phase of 2
                //
                // Let's say we are at position 1. Then the distance is
                // (6-1)-1 = 4. We subtract the phase of 2 which gives us 2.
                // 2 mod 3 is 2.
                // Let's say we are at position 5. Then we calculate the distance
                // (6-1)-5 = 0. We subtract the phase 2 so we are at -2.
                // -2 mod 3 is 1.
                let relative_frame = mod3(
                    cds.range.stop as isize
                                       -1 // because the end position is exclusive
                                       - phase as isize // move to a position that is at a codon start (on the reverse strand)
                                       - genomic_position as isize,
                ); // get the distance and do mod3
                debug_assert!(relative_frame < 3, "relative_frame = {}", relative_frame);

                let mut nucs = String::with_capacity(3);
                //let pos_of_codon_start = middle_index + relative_frame;
                let codon_range = match relative_frame {
                    //                                         v--- because end is eclusive
                    0 => (middle_index - 2..middle_index + 0 + 1),
                    1 => (middle_index - 1..middle_index + 1 + 1),
                    2 => (middle_index - 0..middle_index + 2 + 1),
                    _ => panic!("Modulo arithmethic is broken"),
                };
                for pos in codon_range.rev() {
                    //note the .rev() because we want to treat the final codon "normally"
                    nucs.push(complement(seq[pos]))
                }
                (nucs, relative_frame)
            }
        };

        let new_codon = {
            let mut codon = String::with_capacity(3);
            for (i, c) in old_codon.chars().enumerate() {
                codon.push(if i == relative_frame {
                    match self.gene_info.strand {
                        Strand::Plus => other_nucleotide,
                        Strand::Minus => complement(other_nucleotide),
                    }
                } else {
                    c // same as old codon
                })
            }
            codon
        };
        assert_eq!(old_codon.len(), 3);
        assert_eq!(new_codon.len(), 3);
        let old_aa = translate(&old_codon).unwrap();
        let new_aa = translate(&new_codon).unwrap();

        if old_aa == new_aa {
            MutationType::Synonymous
        } else if old_aa == AminoAcid::Stop {
            MutationType::StopLoss
        } else if new_aa == AminoAcid::Stop {
            MutationType::Nonsense
        } else {
            MutationType::Missense
        }
    }
}

#[derive(Debug)]
pub struct SetupError {
    filename: String,
    message: String,
}

impl fmt::Display for SetupError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}: {}", self.filename, self.message)
    }
}
