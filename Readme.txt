Author: Konstantinos Prousalis
Email: kprousalis@csd.auth.gr

PURPOSE
This program is designed to detect perfect or imperfect diagonal lines for dot matrix
plots. As accepted files are only fasta format.
Only pairwise alignment is allowed for different sequences covering the 
most common patterns:
1. perfect match
2. homologous in the same diagonal
3. homologous that are not in the same diagonal
4. indels by requesting via dotplot
Plotting a sequence by itself may form patterns that are not supported.

HINT 1
The user can download the readily available fasta files from genomic 
databases, save them as querySeq.fasta and refSeq.fasta and run their content.
Beware, remove some initial instructions.

HINT 2
The user can use the seqdotplot(seq1,seq2) function throughout the code
to enable the view of the dotplot per window.

HINT 3
Runtime may suffer an exponential slow down for very large windows 
(better to use windows smaller than 1000) or wide scanning around the 
main diagonal (critical parameter hw, see the code).
