# con-hi

“Con-hi” means “**con**sensus-**hi**ghlighter”.

Latest version is `3.3.c` (2024-12-22 edition).

## Description

This program annotates low-coverage and high-coverage regions of sequences in fasta format using read mapping in BAM format.

### Input

1. Target sequence(s) in fasta format.
2. Read mapping in a **sorted** BAM file.
3. Coverage threshold(s) for searching low-coverage and high-coverage regions.

### Output

- A GenBank file with annotated low-coverage and high-coverage regions.

## Dependencies

- [Python](https://www.python.org/downloads/) 3.6 or later.
- [Biopython](https://biopython.org/) package.
- [samtools](https://github.com/samtools/samtools) 1.13 or later is recommended. Versions from 1.11 to 1.12 are acceptable, but might calculate coverage inaccurately.

You can install Biopython with following command:
```bash
pip3 install biopython
```

You can install samtools by downloading latest release from [samtools page](https://github.com/samtools/samtools) on github. Then follow instrunctions in downloaded INSTALL file.

## Usage

Basic usage is:
```bash
./con-hi.py -f <TARGET_FASTA> -b <MAPPING_BAM>
```

You can specify custom coverage theshold(s) by passing comma-separated list of thresholds with options `-c` and `-C`. For example, following command will annotate:

- regions with coverage below 25 and all regions below 55 (and also with zero coverage);

- regions with coverage greater than 1.5×*M* and greater than 2.0×*M*, where *M* is median coverage.

```
./con-hi.py \
    -f my_sequence.fasta -b my_mapping.sorted.bam \
    -c 25,55 -C 1.5,2.0
```

### Options

```
-f or --target-fasta: *
    File of target sequence(s) in fasta format.

-b or --bam: *
    Sorted BAM file which contains mapping on target sequence(s).

-o or --outfile:
    Output file.
    Deault value: 'highlighted_sequence.gbk'.

-r or --target-seq-ids:
    Comma-separated list of target sequence IDs to process.
    Examples: "seq_1" or "seq_1,seq_9,seq_12".
    Dasta sequence id is the part of its header before the first space.
    Default: process all target sequences.

-c or --lower-coverage-thresholds:
    Comma-separated list of lower coverage threshold(s).
    Default: 10.
    To disable it, specify "-c off", and low-coverage regions won't be annotated.

-n or --no-zero-output:
    Disable annotation of zero-coverage regions.
    Disabled by default.

-C or --upper-coverage-coefficients:
    Comma-separated list of coverage coefficient(s).
    To annotate regions with coverage > 1.7×M,
      where M is median coverage, specify "-C 1.7".
    Default: 2.0.
    To disable it, specify "-C off", and high-coverage regions won't be annotated.

-l or --min-feature-len:
    Minimum length of a feature to output. Must be int >= 0.
    Default: 5 bp.

--circular:
    Target sequence in curcular. Affects only corresponding GenBank field.
    Disabled by default.
    
--organism:
    Organism name. Affects only corresponding GenBank field.
    If it contains spaces, surround it with quotes (see Example 4).
    Empty by default.

-k or --keep-temp-cov-file:
    Don't delete temporary TSV file "coverages.tsv" after work of the program.
    The program creates this file in the same directory where the "-o" file is located.
    Default behaviour is to delete this file afterwards.
```
\* - mandatory option


## Examples

### Example 1. Basic usage

Annotate file `my_sequence.fasta` with default parameters according to mapping from file `my_mapping.sorted.bam`:

```
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam
```

### Example 2. How to use `-c` option

Annotate regions with coverage below 25, fragments with coverages below 50 and regions with zero coverages:

```bash
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam -c 25,50
```

### Example 3. How to use options `-C`, `-n`

Annotate regions with coverage below 25, fragments with coverages below 50. Disable annotation of zero coverage regions:

```
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam -c 25,50 -n
```

### Example 4. How to use options `--circular`, `--organism`

Specify the name of the organism for output file. The sequence is circular:

```bash
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam \
    --circular --organism "Serratia marcescens"
```

### Example 5. How to turn off annotation of low-coverage regions (`-c off`)

Disable annotation of low-coverage regions (`-c off`). Annotate high-coverage regions with coverage above 1.7×*M* and above 2.4×*M*, where *M* is median coverage:

```bash
./con-hi.py -f my_sequence.fasta -b my_mapping.sorted.bam \
    -c off -C 1.7,2.4
```

### Example 6. How to use `-r` option

Target file `my_sequences.fasta` contains the following sequences:

1) a prokaryotic chromosome (sequence id `chr`);

2) one high-copy plasmid (sequence id `plasmid_H1`);

3) two low-copy plasmids (sequence ids `plasmid_L1` and `plasmid_L2`).

One might expect that the more copies a replicon has the higher is its read coverage. Use coverage threshold of 20 for the chromosome, 50 for the high-copy plasmid, and 5 for low-copy plasmids:

```bash
./con-hi.py -f my_sequences.fasta -b my_mapping.sorted.bam \
    -r chr \
    -c 20

./con-hi.py -f my_sequences.fasta -b my_mapping.sorted.bam \
    -r plasmid_H1 \
    -c 50

./con-hi.py -f my_sequences.fasta -b my_mapping.sorted.bam \
    -r plasmid_L1,plasmid_L2 \
    -c 5
```
