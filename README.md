# SCARR

SCARR (v. 1.0) is a Python 2.7 script that takes as input the output from blastn (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) in tabular format ("-outfmt 6" in the command line version of blastn).

From the output blastn file, SCARR tries to identify the best possible combination of alignments to explain the original sequences as DNA rearrangements. It first checks whether the sequence can be explained as a single alignment. If this is not the case, it tests all possible combinations of two alignments, giving each one a score. If the best score passes the threshold specified, the sequence is identified as a rearrangement and SCARR will move on to the next sequence. If no combination of two alignments passes the threshold, SCARR then tests all possible combinations of three alignments and gives them a score, as for the two-alignment search. If still no combination passes the threshold, the sequence is determined unlabelled and SCARR moves on to the next. 

SCARR produces 5 output files:

'AnalysisSummary.txt' is a tab-separated table that contains the result for each sequence identifier found in the blastn output file, along with the best scores for single-alignment, two-alignment, and three-alignment analyses.

'JunctionStats.txt' is a tab-separated table that indicates how many two-alignment (Single Junctions) and three-alignment (Double Junctions) rearrangements are found. Each number is also broken down as rearrangements that have 5 bases or more of homology at their breakpoint junction (Microhomology) or 4 bases or less (No microhomology). Finally, the number of short-range inversions is specified (U-Turns).

'SingleJunctions.txt' is a tab-separated table that lists each sequence identified as a two-alignment rearrangement over 2 lines. The first line contains the blastn information of the alignment with the smallest position on the original sequence, along with the SCARR score and whether the rearrangement contains a microhomology, no microhomology, or a U-turn. The second line contains the blastn information of the alignment with the highest position on the original sequence.

'DoubleJunctions.txt' is a tab-separated table that lists each sequence identified as a three-alignment rearrangement over 3 lines. The first line contains the blastn information of the alignment with the smallest position on the original sequence, along with the SCARR score and whether the first rearrangement junction contains a microhomology, no microhomology, or a U-turn. The second line contains the blastn information of the alignment with the next lowest position on the original sequence, along with the SCARR score and whether the second rearrangement junction contains a microhomology, no microhomology, or a U-turn. The third line contains the blastn information of the alignment with the highest position on the original sequence.

'UnlabelledReads.txt' is a tab-separated table that contains the identifiers of sequences that were not explained by a single alignment, two-alignment rearrangement or three-alignment rearrangement, along with the best SCARR scores for each possibility.

SCARRClassify (v. 1.0) is a Python 2.7 script that takes as input a SingleJunctions.txt file output from SCARR.

SCARRClassify looks at each single junction rearrangement and determines whether it corresponds to a short deletion or duplication (Slippage), the length of homology at the breakpoint junction (Overlap), how many bases were inserted at the breakpoint junction (Insertion), whether the rearrangement is an inversion, deletion, duplication or translocation. If the rearrangement is a deletion or duplication, it measures the length of the deleted or duplicated segment, and provides the position of the first base deleted or duplicated (Breakpoint 1) and the position of the last base deleted or duplicated (Breakpoint 2). For inversions, it measures the distance between the two breakpoints minus any homology, and provides the breakpoint on the + strand (Breakpoint 1) and on the - strand (Breakpoint 2).

SCARRClassify also produces 7 output files:

'SingleJunctionsDetail.txt' is a tab-separated table in which the additional information described in the previous paragraph is added to 'SingleJunctions.txt' 

'Deletions.txt', 'Duplications.txt' and 'Inversions.txt' are tab-separated files that contain information about the rearrangements identified for each type (chromosome, breakpoint positions, length, homology and inserted bases).

'DeletionsGrouped.txt','DuplicationsGrouped.txt' and 'InversionsGrouped.txt' are tab-separated files that contain the same information as 'Deletions.txt', 'Duplications.txt' and 'Inversions.txt', but in which identical rearrangements are merged. Rearrangements are sorted in order of how many times they were identified.
