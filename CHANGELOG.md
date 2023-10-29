# consensus-highlighter changelog

# 2023-10-29 edition

Version `3.1.a`

1. Fix a bug that would cause the program to terminate if `-C` threshold is enabled and a high-coverage region starts at position 0 or ends at position (LENGTH-1) of the reference sequence.

2. Now, the program prints its own warning if length of an output sequence name is too long for "pretty" GenBank representation, according to GenBank standard (Dec 15, 2018) 229.0. In previous versions of the program, BioPython warning used to be emitted, which is not very informative: "Increasing length of locus line to allow long name. This will result in fields that are not in usual positions.".

## 2023-05-25 edition

Version `3.0.b`

Fix a bug preventing `consensus-highlighter` from parsing `samtools version` output correctly if `samtools` is compiled with a flag `-ffile-prefix-map`. In that case, the `samtools version` output contains some non-utf8 characters.

## 2022-07-26 edition

Version `3.0.a`

1. Add option `-l/min-feature-len`. It sets minimum length of an output feature.

2. Add option `-C/upper-coverage-coefficients`. It sets threshold(s) for annotating **high-coverage** regions. For example, to annotate regions with coverage > 1.7×*M*, where *M* is median coverage, you should specify `-C 1.7`. You can specify multiple coefficients: `-C 1.5,1.7`, in the same way as for option `-c`.

3. Change long option name: `-c/coverage-thresholds` -> `-c/lower-coverage-thresholds`.

4. Options `-c` and `-C` can be disabled now: specify `-c off`, `-C off`, and low-coverage or high-coverage regions won't be annotated, respectively.

5. Add option `-k/--keep-temp-cov-file`. If it is specified, temporary file `coverages.tsv` won't be deleted after work of the program.

## 2022-04-26 edition

Version `2.3.a`

1. Add recommendation "samtools `1.13` or later is recommended". This is the version, in which `samtools depth` [had beed](https://github.com/samtools/samtools/releases/tag/1.13) completely rewritten. Since 1.13, `samtools depth` calculates coverage more accurately.

2. Now consensus-highlighter does not crash if samtools version is of the following format: `1.15.1` (three dot-separated numbers). Previously, only two dot-separated numbers were permitted.

## 2022-02-23 edition

Version `2.2.b`

Changes:
1. Now consensus-highlighter removes its temporary file `coverages.tsv`, where coverage value of each base is stored.
2. A bug was fixed that would cause the program to write error message about an unmet dependency to stout instead of stderr.

## 2022-02-10 edition

Version `2.2.a`

Now consensus-highlighter adds a comment to output GenBank files. Here is the example of such a comment:

```
COMMENT     ##Coverage-Data-START##
            Minimum Coverage    :: 0
            Average Coverage    :: 147.54
            Median Coverage     :: 70
            Maximum Coverage    :: 2009
            Zero-coverage bases :: 278 bp
            ##Coverage-Data-END##
```

## 2022-01-24 edition

Version `2.1.a`

Now consensus-highlighter is compatible with samtools 1.13+.

## 2021-12-23 edition

Version `2.0.b`

Now consensus-highlighter calculates and prints average coverage.

## 2021-12-09 edition

Version `2.0.a`

Removed options `-o/--outdir` and `--prefix`.

Added option `-o/--outfile`. And now consensus-highlighter.py writes all GenBank output records to this single output GenBank file.

## 2021-06-15 edition

Version `1.1.a`

`consensus-highlighter` no more piles up coverage features with identical locations. It means that you will not see both "zero coverage" and "coverage < 10" features starting at the same positions and ending at the same positions.

## 2021-06-11 edition

Version `1.0.d`

Added warning messages for following cases:

1. If the program cannot find ids of sequence(s) from `-f` fasta file in the coverage file.
2. If length of sequence in `-f` fasta file is not equal to number of coverage positions reported by `samtools depth` and stored in the coverage file.

## 2021-06-10 edition

Version `1.0.b`

Fixed bug that would cause the program to stumble on a non-extant directory in PATH while checking dependencies.

And then...

Version `1.0.c`

Fixed bug that cause the program to stumble on lowercase input sequences.

## 2021-06-03 edition

Init release. Version `1.0.a`
