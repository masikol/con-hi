# consensus-highlighter changelog

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
