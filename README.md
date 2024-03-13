[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10115908.svg)](https://doi.org/10.5281/zenodo.10115908)

# SWOffinder


SWOffinder is a method based on Smith-Waterman alignment to find all off-target sites up to some edit distance efficiently.
We implemented a trace-back method to find off-target sites under versatile criteria, such as separate limits on the number of insertions, deletions, and mismatches.

# Usage
1. First, you need to compile and build the SmithWatermanOffTarget package: `javac -d bin SmithWatermanOffTarget/*.java`
2. To search, execute the command line `java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign <Genome reference path> <sgRNA list path> <Output path prefix> <maxE> <maxM> <maxMB> <maxB> <Threads> <Best in-a-window> <Best in-a-window size> <PAM> <Allow PAM edits>`

The arguments specification:

1. **Genome reference path**: The path of the FASTA file containing the Genome reference.
2. **sgRNA list**: The path of a text file containing a list of sgRNAs with their PAM (see sgRNAs.txt file for example).
3. **Output path prefix**: The path for the output files with file prefix. files are saved in "CSV" format.
4. **maxE**: Max edits allowed (integer).
5. **maxM**: Max mismatches allowed without bulges (integer).
6. **maxMB**: Max mismatches allowed with bulges (integer).
7. **maxB**: Max bulges allowed (integer).
8. **Threads**: The number of threads to use for the run (integer). 
9. **Best in-a-window**: Flag whether to choose the best off-target site in a window or not (true or false).
10. **Best in-a-window size**: The window size for choosing the best in a window (integer). Please insert even if **Best in-a-window** is false.
11. **PAM**: The PAM sequence (for example, NGG).
12. **Allow PAM edits**: Flag whether to allow PAM edits or not (true or false).


# Requirements:
The code was tested with openjdk version "17.0.3"

