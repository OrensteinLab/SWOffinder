# SWOffinder


SWOffinder, a method based on Smith-Waterman alignment to efficiently find all off-target sites up to some edit distance.
We implemented a trace-back method to find off-target sites under versatile criteria, such as separate limits on the number of insertions, deletions, and mismatches.

# Usage
1. First you need to compile and build the SmithWatermanOffTarget package: `javac -d bin SmithWatermanOffTarget/*.java`
2. To search execute the command line `java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign <Genome reference path> <sgRNA list path> <Output path prefix> <maxE> <maxM> <maxMB> <maxB> <Threads> <Best in-a-window> <Best in-a-window size> <PAM> <Allow PAM edits>`

The arguments specifiction:

1. **Genome reference path**: The path of the FASTA file containing the Genome reference.
2. **sgRNA list**: The path of a text file containing a list of sgRNAs with their PAM (see sgRNAs.txt file for example).
3. **Output path prefix**: The path for the output files with file prefix. files are saved in "CSV" format.
4. **maxE**: Max edits allowed (integer).
5. **maxM**: Max mismatches allowed without bulges (integer).
6. **maxMB**: Max mismatches allowed with bulges (integer).
7. **maxB**: Max bulges allowed (integer).
8. **Threads**: The number of threads to use for the run (integer). 
9. **Best in-a-window**: Flag whether to choose the best off-target site in a window or not (true or false).
10. **Best in-a-window size**: In case of best in-a-window is true, what is the the size of the window (integer).
11. **PAM**: The PAM seuqnce (for example, NGG).
12. **Allow PAM edits**: Flag whether to allow PAM edits or no (true of false).


# Requirements:
The code was tested with openjdk version "17.0.3"

