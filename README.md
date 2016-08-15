# Execute FMextend error correction
FMextend is a new algorithm for error correction in StriDe correct step.

Command: 

stride correct -a fmextend -t 30 -k 31 -K 15 -x 3 reads.fa -o READ.ECOLr.fasta

-a : the corretion algorithm.

-k : the kmer size for checking kmer solid or not.

-K : the kmer size when align in extend step.
