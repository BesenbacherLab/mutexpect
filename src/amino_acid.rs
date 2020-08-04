#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AminoAcid {
    Alanine,
    Arginine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    GlutamicAcid,
    Glutamine,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    Stop,
}

const GENETIC_CODE: [AminoAcid; 64] = [
    AminoAcid::Lysine,        //AAA
    AminoAcid::Asparagine,    //AAC
    AminoAcid::Lysine,        //AAG
    AminoAcid::Asparagine,    //AAT
    AminoAcid::Threonine,     //ACA
    AminoAcid::Threonine,     //ACC
    AminoAcid::Threonine,     //ACG
    AminoAcid::Threonine,     //ACT
    AminoAcid::Arginine,      //AGA
    AminoAcid::Serine,        //AGC
    AminoAcid::Arginine,      //AGG
    AminoAcid::Serine,        //AGT
    AminoAcid::Isoleucine,    //ATA
    AminoAcid::Isoleucine,    //ATC
    AminoAcid::Methionine,    //ATG
    AminoAcid::Isoleucine,    //ATT
    AminoAcid::Glutamine,     //CAA
    AminoAcid::Histidine,     //CAC
    AminoAcid::Glutamine,     //CAG
    AminoAcid::Histidine,     //CAT
    AminoAcid::Proline,       //CCA
    AminoAcid::Proline,       //CCC
    AminoAcid::Proline,       //CCG
    AminoAcid::Proline,       //CCT
    AminoAcid::Arginine,      //CGA
    AminoAcid::Arginine,      //CGC
    AminoAcid::Arginine,      //CGG
    AminoAcid::Arginine,      //CGT
    AminoAcid::Leucine,       //CTA
    AminoAcid::Leucine,       //CTC
    AminoAcid::Leucine,       //CTG
    AminoAcid::Leucine,       //CTT
    AminoAcid::GlutamicAcid,  //GAA
    AminoAcid::AsparticAcid,  //GAC
    AminoAcid::GlutamicAcid,  //GAG
    AminoAcid::AsparticAcid,  //GAT
    AminoAcid::Alanine,       //GCA
    AminoAcid::Alanine,       //GCC
    AminoAcid::Alanine,       //GCG
    AminoAcid::Alanine,       //GCT
    AminoAcid::Glycine,       //GGA
    AminoAcid::Glycine,       //GGC
    AminoAcid::Glycine,       //GGG
    AminoAcid::Glycine,       //GGT
    AminoAcid::Valine,        //GTA
    AminoAcid::Valine,        //GTC
    AminoAcid::Valine,        //GTG
    AminoAcid::Valine,        //GTT
    AminoAcid::Stop,          //TAA
    AminoAcid::Tyrosine,      //TAC
    AminoAcid::Stop,          //TAG
    AminoAcid::Tyrosine,      //TAT
    AminoAcid::Serine,        //TCA
    AminoAcid::Serine,        //TCC
    AminoAcid::Serine,        //TCG
    AminoAcid::Serine,        //TCT
    AminoAcid::Stop,          //TGA
    AminoAcid::Cysteine,      //TGC
    AminoAcid::Tryptophan,    //TGG
    AminoAcid::Cysteine,      //TGT
    AminoAcid::Leucine,       //TTA
    AminoAcid::Phenylalanine, //TTC
    AminoAcid::Leucine,       //TTG
    AminoAcid::Phenylalanine, //TTG
];

pub fn translate(codon: &str) -> Option<AminoAcid> {
    if codon.len() == 3 {
        let mut index = 0;
        for c in codon.chars() {
            let v = match c {
                // translate to base 4
                'A' | 'a' => 0,
                'C' | 'c' => 1,
                'G' | 'g' => 2,
                'T' | 't' | 'U' | 'u' => 3,
                _ => return None,
            };
            index <<= 2; // times 4
            index += v;
        }
        Some(GENETIC_CODE[index])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_translate() {
        assert_eq!(translate(""), None);
        assert_eq!(translate("A"), None);
        assert_eq!(translate("AT"), None);
        assert_eq!(translate("UAGA"), None);
        assert_eq!(translate("NOPE"), None);

        //full genetic code test
        assert_eq!(translate("AAA"), Some(AminoAcid::Lysine));
        assert_eq!(translate("AAC"), Some(AminoAcid::Asparagine));
        assert_eq!(translate("AAG"), Some(AminoAcid::Lysine));
        assert_eq!(translate("AAT"), Some(AminoAcid::Asparagine));
        assert_eq!(translate("ACA"), Some(AminoAcid::Threonine));
        assert_eq!(translate("ACC"), Some(AminoAcid::Threonine));
        assert_eq!(translate("ACG"), Some(AminoAcid::Threonine));
        assert_eq!(translate("ACT"), Some(AminoAcid::Threonine));
        assert_eq!(translate("AGA"), Some(AminoAcid::Arginine));
        assert_eq!(translate("AGC"), Some(AminoAcid::Serine));
        assert_eq!(translate("AGG"), Some(AminoAcid::Arginine));
        assert_eq!(translate("AGT"), Some(AminoAcid::Serine));
        assert_eq!(translate("ATA"), Some(AminoAcid::Isoleucine));
        assert_eq!(translate("ATC"), Some(AminoAcid::Isoleucine));
        assert_eq!(translate("ATG"), Some(AminoAcid::Methionine));
        assert_eq!(translate("ATT"), Some(AminoAcid::Isoleucine));
        assert_eq!(translate("CAA"), Some(AminoAcid::Glutamine));
        assert_eq!(translate("CAC"), Some(AminoAcid::Histidine));
        assert_eq!(translate("CAG"), Some(AminoAcid::Glutamine));
        assert_eq!(translate("CAT"), Some(AminoAcid::Histidine));
        assert_eq!(translate("CCA"), Some(AminoAcid::Proline));
        assert_eq!(translate("CCC"), Some(AminoAcid::Proline));
        assert_eq!(translate("CCG"), Some(AminoAcid::Proline));
        assert_eq!(translate("CCT"), Some(AminoAcid::Proline));
        assert_eq!(translate("CGA"), Some(AminoAcid::Arginine));
        assert_eq!(translate("CGC"), Some(AminoAcid::Arginine));
        assert_eq!(translate("CGG"), Some(AminoAcid::Arginine));
        assert_eq!(translate("CGT"), Some(AminoAcid::Arginine));
        assert_eq!(translate("CTA"), Some(AminoAcid::Leucine));
        assert_eq!(translate("CTC"), Some(AminoAcid::Leucine));
        assert_eq!(translate("CTG"), Some(AminoAcid::Leucine));
        assert_eq!(translate("CTT"), Some(AminoAcid::Leucine));
        assert_eq!(translate("GAA"), Some(AminoAcid::GlutamicAcid));
        assert_eq!(translate("GAC"), Some(AminoAcid::AsparticAcid));
        assert_eq!(translate("GAG"), Some(AminoAcid::GlutamicAcid));
        assert_eq!(translate("GAT"), Some(AminoAcid::AsparticAcid));
        assert_eq!(translate("GCA"), Some(AminoAcid::Alanine));
        assert_eq!(translate("GCC"), Some(AminoAcid::Alanine));
        assert_eq!(translate("GCG"), Some(AminoAcid::Alanine));
        assert_eq!(translate("GCT"), Some(AminoAcid::Alanine));
        assert_eq!(translate("GGA"), Some(AminoAcid::Glycine));
        assert_eq!(translate("GGC"), Some(AminoAcid::Glycine));
        assert_eq!(translate("GGG"), Some(AminoAcid::Glycine));
        assert_eq!(translate("GGT"), Some(AminoAcid::Glycine));
        assert_eq!(translate("GTA"), Some(AminoAcid::Valine));
        assert_eq!(translate("GTC"), Some(AminoAcid::Valine));
        assert_eq!(translate("GTG"), Some(AminoAcid::Valine));
        assert_eq!(translate("GTT"), Some(AminoAcid::Valine));
        assert_eq!(translate("TAA"), Some(AminoAcid::Stop));
        assert_eq!(translate("TAC"), Some(AminoAcid::Tyrosine));
        assert_eq!(translate("TAG"), Some(AminoAcid::Stop));
        assert_eq!(translate("TAT"), Some(AminoAcid::Tyrosine));
        assert_eq!(translate("TCA"), Some(AminoAcid::Serine));
        assert_eq!(translate("TCC"), Some(AminoAcid::Serine));
        assert_eq!(translate("TCG"), Some(AminoAcid::Serine));
        assert_eq!(translate("TCT"), Some(AminoAcid::Serine));
        assert_eq!(translate("TGA"), Some(AminoAcid::Stop));
        assert_eq!(translate("TGC"), Some(AminoAcid::Cysteine));
        assert_eq!(translate("TGG"), Some(AminoAcid::Tryptophan));
        assert_eq!(translate("TGT"), Some(AminoAcid::Cysteine));
        assert_eq!(translate("TTA"), Some(AminoAcid::Leucine));
        assert_eq!(translate("TTC"), Some(AminoAcid::Phenylalanine));
        assert_eq!(translate("TTG"), Some(AminoAcid::Leucine));
        assert_eq!(translate("TTT"), Some(AminoAcid::Phenylalanine));

        // test some RNA codons
        assert_eq!(translate("AUU"), Some(AminoAcid::Isoleucine));
        assert_eq!(translate("GCU"), Some(AminoAcid::Alanine));
        assert_eq!(translate("UUU"), Some(AminoAcid::Phenylalanine));
        assert_eq!(translate("UAU"), Some(AminoAcid::Tyrosine));
    }
}
