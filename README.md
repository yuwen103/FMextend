# Execute FMextend error correction
FMextend is a new algorithm for error correction step in StriDe.

StriDe: https://github.com/ythuang0522/StriDe

Command: 

stride correct -a fmextend -t 30 -k 31 -K 15 InputName.fa -o OutputName.fa

InputName.fa : Input file(fasta)

OutputName.fa : the file of corrected read.

-a : the corretion algorithm.

-k : the kmer size for checking kmer solid or not.

-K : the kmer size when align in extend step.
