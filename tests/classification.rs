//! Classification of possible mutations by `for_each_possible_mutation`.
//!
//! Hand-checked cases on small genes, a comparison with an independent
//! reference (spliced codons translated with the standard genetic code) on
//! random plus-strand genes, and a symmetry test: a minus-strand gene must be
//! classified exactly like its reverse complement on the plus strand.
//!
//! Conventions (as in `SeqAnnotation`): 0-based half-open intervals; the CDS
//! phase is the number of bases at the CDS's 5' end (on the coding strand)
//! before the first complete codon; codons split between two CDSs are
//! classified as `Unknown`.

use std::collections::BTreeMap;

use mutexpect::interval::Interval;
use mutexpect::{
    for_each_possible_mutation, Alteration, MutationType, Phase, SeqAnnotation, Strand, CDS,
};

use MutationType::*;

/// (position, alteration) -> type; alterations: a base, 'I' (in-frame indel)
/// or 'F' (frameshift indel).
type Classes = BTreeMap<(usize, char), MutationType>;

fn phase(p: u8) -> Phase {
    match p {
        0 => Phase::Zero,
        1 => Phase::One,
        _ => Phase::Two,
    }
}

fn annotation(
    strand: Strand,
    range: (usize, usize),
    exons: &[(usize, usize)],
    cds: &[(usize, usize, u8)],
) -> SeqAnnotation {
    SeqAnnotation::new(
        "T1".to_string(),
        "chrT".to_string(),
        Interval::new(range.0, range.1).unwrap(),
        strand,
        exons
            .iter()
            .map(|&(s, e)| Interval::new(s, e).unwrap())
            .collect(),
        cds.iter()
            .map(|&(s, e, p)| CDS::new(Interval::new(s, e).unwrap(), phase(p)))
            .collect(),
    )
}

/// Classify every possible mutation of `ann` in `genome` (all types kept).
fn classify(genome: &str, ann: &SeqAnnotation, flank: usize, filter_plof: bool) -> Classes {
    // Like genovo: one base more on the right (the walk also visits range.stop).
    let seq = &genome[ann.range.start - flank..ann.range.stop + flank + 1];
    let mut out = Classes::new();
    for_each_possible_mutation(seq, ann, flank, true, true, filter_plof, |site| {
        assert_eq!(
            site.ref_base.name(),
            genome.as_bytes()[site.genomic_position] as char
        );
        for m in site.mutations {
            let alt = match m.alteration {
                Alteration::Point(b) => b.name(),
                Alteration::InFrameIndel => 'I',
                Alteration::FrameshiftIndel => 'F',
            };
            assert!(out.insert((site.genomic_position, alt), m.mutation_type).is_none());
        }
        Ok::<(), ()>(())
    })
    .unwrap();
    out
}

fn get(c: &Classes, pos: usize, alt: char) -> Option<MutationType> {
    c.get(&(pos, alt)).copied()
}

// --- hand-checked cases ------------------------------------------------------

/// Single-exon plus-strand gene: ATG TGG TAC TAA at positions 5..17.
#[test]
fn plus_strand_codons() {
    //            0    5  8  11 14 17
    let genome = "GGGGGATGTGGTACTAAGGGGG";
    let ann = annotation(Strand::Plus, (5, 17), &[(5, 17)], &[(5, 17, 0)]);
    let c = classify(genome, &ann, 2, false);
    // TGG (Trp): third base G -> A gives TGA (stop), -> T/C gives Cys.
    assert_eq!(get(&c, 10, 'A'), Some(Nonsense));
    assert_eq!(get(&c, 10, 'T'), Some(Missense));
    assert_eq!(get(&c, 10, 'C'), Some(Missense));
    // TAC (Tyr): C -> T is TAT (Tyr), -> A/G is TAA/TAG (stop).
    assert_eq!(get(&c, 13, 'T'), Some(Synonymous));
    assert_eq!(get(&c, 13, 'A'), Some(Nonsense));
    assert_eq!(get(&c, 13, 'G'), Some(Nonsense));
    // TAA (stop): second base A -> G is TGA (still stop), -> C is TCA (Ser).
    assert_eq!(get(&c, 15, 'G'), Some(Synonymous));
    assert_eq!(get(&c, 15, 'C'), Some(StopLoss));
    // Any change of the ATG start codon is a loss of the start codon.
    for p in 5..8 {
        for alt in ['A', 'C', 'G', 'T'] {
            if alt != genome.as_bytes()[p] as char {
                assert_eq!(get(&c, p, alt), Some(StartCodon), "{p} {alt}");
            }
        }
    }
    // Indels are anchored at CDS positions after the first one.
    assert_eq!(get(&c, 5, 'I'), None);
    for p in 6..17 {
        assert_eq!(get(&c, p, 'I'), Some(InFrameIndel), "{p}");
        assert_eq!(get(&c, p, 'F'), Some(FrameshiftIndel), "{p}");
    }
    // Outside the CDS (flanks are not part of the region): nothing.
    assert!(c.keys().all(|&(p, _)| (5..18).contains(&p)));
    // Position 17 (one past the region) is not in an exon: Unknown.
    assert_eq!(get(&c, 17, 'A'), Some(Unknown));
}

#[test]
fn start_codon_loss_needs_an_atg_in_phase_0() {
    // First codon is not ATG (e.g. an incomplete 5' annotation): ordinary
    // codon classification. CTG -> CTA (Leu) is synonymous.
    let genome = "GGGGGCTGTGGTACTAAGGGGG";
    let ann = annotation(Strand::Plus, (5, 17), &[(5, 17)], &[(5, 17, 0)]);
    let c = classify(genome, &ann, 2, false);
    assert_eq!(get(&c, 7, 'A'), Some(Synonymous));
    assert_eq!(get(&c, 5, 'A'), Some(Missense)); // ATG, but not a loss of a start codon
    // ATG at the CDS start, but the first CDS has phase 1: no start codon.
    let genome = "GGGGGATGTGGTACTAAGGGGG";
    let ann = annotation(Strand::Plus, (5, 17), &[(5, 17)], &[(5, 17, 1)]);
    let c = classify(genome, &ann, 2, false);
    assert!(c.values().all(|&t| t != StartCodon));
    // Minus strand: the start codon is at the high end of the last CDS
    // (CAT on the plus strand at positions 14..17).
    let genome = "CCCCCTTAGTACCACATCCCCC";
    let ann = annotation(Strand::Minus, (5, 17), &[(5, 17)], &[(5, 17, 0)]);
    let c = classify(genome, &ann, 2, false);
    for p in 14..17 {
        assert_eq!(get(&c, p, 'G').or(get(&c, p, 'C')), Some(StartCodon), "{p}");
    }
    let point = |a: char| !matches!(a, 'I' | 'F');
    assert!(c
        .iter()
        .filter(|((_, a), _)| point(*a))
        .all(|(&(p, _), &t)| (t == StartCodon) == (14..17).contains(&p)));
}

#[test]
fn filter_plof_drops_nonsense_and_frameshifts_in_last_50_coding_bases() {
    // 20 codons GCT (Ala) with a TGG at codon 2 and codon 19: the CDS is 60
    // bases, so the last 50 coding bases start at position 5 + 10 = 15.
    let mut cds = "GCT".repeat(20);
    cds.replace_range(3..6, "TGG");
    cds.replace_range(54..57, "TGG");
    let genome = format!("GGGGG{cds}GGGGG");
    let ann = annotation(Strand::Plus, (5, 65), &[(5, 65)], &[(5, 65, 0)]);
    let c = classify(&genome, &ann, 2, false);
    let f = classify(&genome, &ann, 2, true);
    // TGG -> TGA at codon 2 (position 10) stays nonsense; at codon 19
    // (position 61) it is in the last 50 coding bases.
    assert_eq!(get(&c, 10, 'A'), Some(Nonsense));
    assert_eq!(get(&f, 10, 'A'), Some(Nonsense));
    assert_eq!(get(&c, 61, 'A'), Some(Nonsense));
    assert_eq!(get(&f, 61, 'A'), Some(Unknown));
    // Frameshifts: kept before position 15, dropped from 15 on; in-frame
    // indels are always kept.
    assert_eq!(get(&f, 14, 'F'), Some(FrameshiftIndel));
    assert_eq!(get(&f, 15, 'F'), None);
    assert_eq!(get(&f, 15, 'I'), Some(InFrameIndel));
    // Nothing else changes.
    let changed: Vec<_> = c.iter().filter(|(k, v)| f.get(k) != Some(v)).collect();
    assert_eq!(changed.len(), 2 + (64 - 15 + 1), "{changed:?}");
}

/// Two exons with a canonical intron GT...AG and a codon split between them.
#[test]
fn splice_sites_and_split_codons() {
    // exon 1: 5..12 = ATG GCT G   (phase 0; the last base starts a split codon)
    // intron: 12..20 = GTAAAAAG
    // exon 2: 20..30 = CT TGG TAA C (phase 2: CT completes the split codon)
    let genome = "GGGGGATGGCTGGTAAAAAGCTTGGTAACGGGGG";
    let ann = annotation(
        Strand::Plus,
        (5, 30),
        &[(5, 12), (20, 30)],
        &[(5, 12, 0), (20, 29, 2)],
    );
    let c = classify(genome, &ann, 2, false);
    for alt in ['A', 'C', 'T'] {
        assert_eq!(get(&c, 12, alt), Some(SpliceSite)); // [G]T
        assert_eq!(get(&c, 19, alt), Some(SpliceSite)); // A[G]
    }
    assert_eq!(get(&c, 13, 'C'), Some(SpliceSite)); // G[T]
    assert_eq!(get(&c, 18, 'C'), Some(SpliceSite)); // [A]G
    for p in 14..18 {
        assert_eq!(get(&c, p, 'C'), Some(Intronic), "{p}");
    }
    // The split codon G|CT is ignored.
    assert_eq!(get(&c, 11, 'A'), Some(Unknown));
    assert_eq!(get(&c, 20, 'A'), Some(Unknown));
    assert_eq!(get(&c, 21, 'A'), Some(Unknown));
    // TGG in exon 2 (positions 22..25) is read in phase.
    assert_eq!(get(&c, 24, 'A'), Some(Nonsense));
    assert_eq!(get(&c, 24, 'T'), Some(Missense));
    // GCT (Ala) in exon 1: third base is four-fold degenerate.
    for alt in ['A', 'C', 'G'] {
        assert_eq!(get(&c, 10, alt), Some(Synonymous));
    }
    // Position 29 is exonic but outside the CDS (UTR).
    assert_eq!(get(&c, 29, 'A'), Some(Unknown));
    // No indels in the intron.
    assert!((12..20).all(|p| get(&c, p, 'I').is_none()));
}

#[test]
fn non_canonical_intron_has_no_splice_sites() {
    // Same as above, but the intron starts CT and ends AC (the minus-strand
    // motifs), which are not splice sites on the plus strand.
    let genome = "GGGGGATGGCTGCTAAAAACCTTGGTAACGGGGG";
    let ann = annotation(
        Strand::Plus,
        (5, 30),
        &[(5, 12), (20, 30)],
        &[(5, 12, 0), (20, 29, 2)],
    );
    let c = classify(genome, &ann, 2, false);
    for p in 12..20 {
        let alt = if &genome[p..p + 1] == "A" { 'C' } else { 'A' };
        assert_eq!(get(&c, p, alt), Some(Intronic), "{p}");
    }
}

#[test]
fn intronic_and_unknown_filters() {
    let genome = "GGGGGATGGCTGGTAAAAAGCTTGGTAACGGGGG";
    let ann = annotation(
        Strand::Plus,
        (5, 30),
        &[(5, 12), (20, 30)],
        &[(5, 12, 0), (20, 29, 2)],
    );
    let seq = &genome[3..33];
    let mut types = Vec::new();
    for_each_possible_mutation(seq, &ann, 2, false, false, false, |site| {
        types.extend(site.mutations.iter().map(|m| m.mutation_type));
        Ok::<(), ()>(())
    })
    .unwrap();
    assert!(!types.contains(&Intronic));
    assert!(!types.contains(&Unknown));
    assert!(types.contains(&SpliceSite));
}

#[test]
fn minus_strand_codons() {
    // Reverse complement of plus_strand_codons: the gene ATG TGG TAC TAA
    // reads TTA GTA CCA CAT on the plus strand at positions 5..17.
    let genome = "CCCCCTTAGTACCACATCCCCC";
    let ann = annotation(Strand::Minus, (5, 17), &[(5, 17)], &[(5, 17, 0)]);
    let c = classify(genome, &ann, 2, false);
    // Third base of TGG (coding strand G) is plus-strand C at position 11.
    assert_eq!(get(&c, 11, 'T'), Some(Nonsense)); // G -> A on the coding strand
    assert_eq!(get(&c, 11, 'A'), Some(Missense));
    // TAA stop: second base is plus-strand T at position 6; T -> C is A -> G.
    assert_eq!(get(&c, 6, 'C'), Some(Synonymous));
    assert_eq!(get(&c, 6, 'G'), Some(StopLoss));
}

// --- reference implementation and random genes -------------------------------

const CODE: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

fn translate(codon: &[u8]) -> u8 {
    let idx = |b: u8| match b {
        b'T' => 0,
        b'C' => 1,
        b'A' => 2,
        b'G' => 3,
        _ => panic!("base {}", b as char),
    };
    CODE[idx(codon[0]) * 16 + idx(codon[1]) * 4 + idx(codon[2])]
}

struct Gene {
    genome: Vec<u8>,
    range: (usize, usize),
    exons: Vec<(usize, usize)>,
    cds: Vec<(usize, usize, u8)>,
}

/// Classes of all point mutations of a plus-strand gene, computed directly
/// from the spliced coding sequence.
fn reference(g: &Gene, filter_plof: bool) -> Classes {
    let mut coding: Vec<usize> = g.cds.iter().flat_map(|&(s, e, _)| s..e).collect();
    coding.sort();
    let last50: Vec<usize> = coding.iter().rev().take(50).copied().collect();
    let mut out = Classes::new();
    for p in g.range.0..=g.range.1 {
        let r = g.genome[p];
        for alt in [b'A', b'C', b'G', b'T'] {
            if alt == r {
                continue;
            }
            let (s0, e0, ph0) = g.cds[0];
            let start_codon =
                ph0 == 0 && e0 - s0 >= 3 && p >= s0 && p < s0 + 3 && &g.genome[s0..s0 + 3] == b"ATG";
            let t = if start_codon {
                StartCodon
            } else if let Some(&(s, e, ph)) = g.cds.iter().find(|c| c.0 <= p && p < c.1) {
                let first = s + ph as usize;
                if p < first {
                    Unknown
                } else {
                    let start = p - (p - first) % 3;
                    if start + 3 > e {
                        Unknown
                    } else {
                        let old = &g.genome[start..start + 3];
                        let mut new = old.to_vec();
                        new[p - start] = alt;
                        let (a, b) = (translate(old), translate(&new));
                        if a == b {
                            Synonymous
                        } else if a == b'*' {
                            StopLoss
                        } else if b == b'*' {
                            if filter_plof && last50.contains(&p) {
                                Unknown
                            } else {
                                Nonsense
                            }
                        } else {
                            Missense
                        }
                    }
                }
            } else if let Some(w) = g.exons.windows(2).find(|w| w[0].1 <= p && p < w[1].0) {
                let (is, ie) = (w[0].1, w[1].0);
                let donor = &g.genome[is..is + 2] == b"GT";
                let acceptor = &g.genome[ie - 2..ie] == b"AG";
                if (donor && (p == is || p == is + 1)) || (acceptor && (p == ie - 2 || p == ie - 1))
                {
                    SpliceSite
                } else {
                    Intronic
                }
            } else {
                Unknown
            };
            out.insert((p, alt as char), t);
        }
    }
    out
}

/// Small deterministic random number generator (no dependencies).
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0 >> 33
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// A random plus-strand gene with 1-4 exons, a CDS spanning a random part of
/// them, consistent phases, and introns that are canonical most of the time.
fn random_gene(rng: &mut Lcg) -> Gene {
    let flank = 6;
    let mut exons = Vec::new();
    let mut pos = flank + 3 + rng.below(5);
    let start = pos;
    let n_exons = 1 + rng.below(4);
    for i in 0..n_exons {
        let len = 4 + rng.below(40);
        exons.push((pos, pos + len));
        pos += len;
        if i + 1 < n_exons {
            pos += 6 + rng.below(10);
        }
    }
    let stop = pos;
    let len = stop + flank + 3;
    let mut genome: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.below(4)]).collect();
    for w in exons.windows(2) {
        if rng.below(4) != 0 {
            genome[w[0].1] = b'G';
            genome[w[0].1 + 1] = b'T';
        }
        if rng.below(4) != 0 {
            genome[w[1].0 - 2] = b'A';
            genome[w[1].0 - 1] = b'G';
        }
    }
    // CDS: from a position in the first exon to a position in the last.
    let cds_start = exons[0].0 + rng.below(3);
    let last = *exons.last().unwrap();
    let cds_stop = last.1 - rng.below(3);
    let mut cds = Vec::new();
    let mut carried = rng.below(3); // phase of the first CDS
    for &(s, e) in &exons {
        let (s, e) = (s.max(cds_start), e.min(cds_stop));
        if s >= e {
            continue;
        }
        let ph = carried.min(e - s);
        cds.push((s, e, ph as u8));
        // bases left over at the 3' end start a codon continued in the next CDS
        let rest = (e - s - ph) % 3;
        carried = (3 - rest) % 3;
    }
    if cds.is_empty() || cds[0].1 - cds[0].0 < 3 {
        return random_gene(rng); // keep the start codon inside the first CDS
    }
    if rng.below(2) == 0 {
        let s = cds[0].0;
        genome[s..s + 3].copy_from_slice(b"ATG");
    }
    Gene {
        genome,
        range: (start, stop),
        exons,
        cds,
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b,
        })
        .collect()
}

fn complement(c: char) -> char {
    match c {
        'A' => 'T',
        'C' => 'G',
        'G' => 'C',
        'T' => 'A',
        x => x,
    }
}

#[test]
fn random_plus_strand_genes_match_reference() {
    let mut rng = Lcg(1);
    for case in 0..300 {
        let g = random_gene(&mut rng);
        let genome = String::from_utf8(g.genome.clone()).unwrap();
        let ann = annotation(Strand::Plus, g.range, &g.exons, &g.cds);
        for plof in [false, true] {
            let ours: Classes = classify(&genome, &ann, 2, plof)
                .into_iter()
                .filter(|((_, a), _)| !matches!(a, 'I' | 'F'))
                .collect();
            let want = reference(&g, plof);
            assert_eq!(ours, want, "case {case}, plof {plof}: {genome}");
        }
    }
}

/// The minus-strand gene obtained by reverse-complementing the genome must be
/// classified like the plus-strand gene (point mutations; indels are anchored
/// on the left in genomic coordinates, so they are not symmetric).
#[test]
fn minus_strand_mirrors_plus_strand() {
    let mut rng = Lcg(2);
    for case in 0..300 {
        let g = random_gene(&mut rng);
        let n = g.genome.len();
        let flip = |(s, e): (usize, usize)| (n - e, n - s);
        let plus = annotation(Strand::Plus, g.range, &g.exons, &g.cds);
        let mut exons: Vec<_> = g.exons.iter().map(|&x| flip(x)).collect();
        exons.reverse();
        let mut cds: Vec<_> = g
            .cds
            .iter()
            .map(|&(s, e, p)| {
                let (s, e) = flip((s, e));
                (s, e, p)
            })
            .collect();
        cds.reverse();
        let minus = annotation(Strand::Minus, flip(g.range), &exons, &cds);
        let gp = String::from_utf8(g.genome.clone()).unwrap();
        let gm = String::from_utf8(revcomp(&g.genome)).unwrap();
        for plof in [false, true] {
            let p = classify(&gp, &plus, 2, plof);
            let m = classify(&gm, &minus, 2, plof);
            for (&(pos, alt), &t) in &p {
                if matches!(alt, 'I' | 'F') || pos == g.range.1 {
                    continue; // indels; the extra position past the region end
                }
                let mirrored = m.get(&(n - 1 - pos, complement(alt))).copied();
                assert_eq!(mirrored, Some(t), "case {case}, plof {plof}, {pos} {alt}");
            }
        }
    }
}

