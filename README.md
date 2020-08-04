# Utility library for determining the probability of point mutations

This library is faily specific for my own purposes. It's on crates.io so that I
can pull it as a dependency for my project. The code is mature at this point.
But it could really use some more doc comments.

This library reads genomic sequences of transcripts and mutates every
nucleotide in that transcript to determine how many mutations of a type
(synonymous, missense, nonsense, ...) can be generated and what their
probabilities are. This is useful to when you try to identify disease genes by
comparing mutated genes across affected individuals and you need to determine
how many mutations you can expect by chance.
