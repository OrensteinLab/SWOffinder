[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10115908.svg)](https://doi.org/10.5281/zenodo.10115908)

# SWOffinder

SWOffinder is a method based on Smith-Waterman alignment to find all off-target sites up to some edit distance efficiently.
We implemented a trace-back method to find off-target sites under versatile criteria, such as separate limits on the number of insertions, deletions, and mismatches.

If you use SWOffinder in your research, please cite our paper:
> [SWOffinder: A Fast and Versatile Tool for CRISPR Off-Target Site Detection](https://www.sciencedirect.com/science/article/pii/S2589004223026342)

# Usage

### 1. Compile
```bash
javac -d bin SmithWatermanOffTarget/*.java
```

### 2. Run
```bash
java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign [options]
```

For a full list of options run:
```bash
java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign --help
```

### Required arguments

| Argument | Description |
|----------|-------------|
| `--genome` | Path to the genome reference FASTA file |
| `--sgrnas` | Path to a text file containing sgRNA sequences with PAM (see `sgRNAs.txt` for example) |
| `--output` | Output path prefix — CSV files will be created here (e.g. `results/my_run_`) |
| `--pam` | PAM sequence (e.g. `NGG`, `NRG`). IUPAC codes are supported |
| `--maxE` | Max edit distance allowed (integer) |
| `--maxM` | Max mismatches allowed without bulges (integer) |
| `--maxMB` | Max mismatches allowed with bulges (integer) |
| `--maxB` | Max bulges allowed (integer) |

### Optional arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--threads` | `8` | Number of threads to use |
| `--bestWindow` | `false` | If `true`, only the best off-target site within a window is reported |
| `--windowSize` | `50` | Window size used when `--bestWindow` is `true` |
| `--allowPamEdits` | `false` | If `true`, mismatches in the PAM region are allowed |
| `--dpWindowSize` | `10000` | Smith-Waterman DP chunk size in base pairs (see performance tips below) |

### Example
```bash
java -cp bin SmithWatermanOffTarget.SmithWatermanOffTargetSearchAlign \
  --genome hg38.fa \
  --sgrnas sgRNAs.txt \
  --output results/my_run_ \
  --pam NGG \
  --maxE 6 \
  --maxM 6 \
  --maxMB 4 \
  --maxB 1
```

> **Note:** sgRNA sequences must include the PAM region (e.g. `GTCCCTAGTGGCCCCACTGTNGG`). A warning will be printed if the PAM is not detected at the end of a sequence.

# Performance Tips

SWOffinder scans a genome-scale sequence, so runtime can be significant. Two parameters have the largest impact on performance:

**`--threads`** — Each thread processes a separate chunk of the genome in parallel. The default is 8, but if your machine has more CPU cores available you can increase this significantly. To check how many cores you have:
- Linux/macOS: `nproc` or `sysctl -n hw.logicalcpu`
- A good rule of thumb is to set `--threads` to the number of logical CPU cores, e.g. `--threads 16` or `--threads 32`.

**`--dpWindowSize`** — Controls the size of each genome chunk processed in a single Smith-Waterman pass. Larger chunks reduce the overhead of thread scheduling and I/O. If your machine has sufficient RAM, increasing this can meaningfully speed up the run. Each thread holds a DP matrix of size `targetLength × dpWindowSize` integers in memory, so the total memory usage is approximately:

```
threads × targetLength × dpWindowSize × 4 bytes
```

For example, with 16 threads, a 23 bp target, and `--dpWindowSize 100000`:
```
16 × 23 × 100000 × 4 bytes ≈ 147 MB
```

A reasonable starting point if you have 16 GB+ RAM is `--dpWindowSize 100000`.

> **Note:** The optimal balance between `--threads` and `--dpWindowSize` depends on your specific machine architecture, memory bandwidth, and CPU cache behavior. We recommend experimenting with different combinations on your machine to find the best performance for your setup.

# Requirements

The code was tested with openjdk version "17.0.3".