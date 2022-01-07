*rhea* conducts SARS-CoV-2 consensus assembly quality checks and if necessary masking.

# rhea
The standard practice for assembling SARS-CoV-2 genomes is to align the sequenced reads to a reference genome ([MN908947.3](https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3)) to create a consensus assembly. In most cases a satisfactory consensus assembly is created, but its often difficult to determine what thresholds were used in the creation of the consensus. *rhea* allows you standardize your consensus assemblies across multiple assemblers. With *rhea* you can set your own coverage cutoffs for sequence masking and identify instances where the reference base was included when it shouldn't have been.

## Installation
*rhea* will be available from Bioconda upon initial release

Temporary installation:
```
# Conda installation
mamba create -n rhea -c conda-forge -c bioconda biopython pysam 'python>3.6' pyvcf
```

# Example Usage
In order to use *rhea* you must provide a few inputs:

1. The sample name used for naming outputs
2. Consensus assembly in FASTA format
3. VCF file containing any mutations observed
4. The BAM file containing read alignments to the make consensus calls

### Usage
```
usage: rhea [-h] [--min_cov INT] [--min_sgene_cov INT] [--width INT] [--outdir STR] [--debug] [--version] SAMPLE FASTA VCF BAM

rhea (v0.0.1) - SARS-CoV-2 consensus assembly quality checks and corrections

positional arguments:
  SAMPLE               Sample name to use for naming outputs
  FASTA                Consensus assembly in FASTA format
  VCF                  VCF file from the consensus assembly
  BAM                  BAM file used to create consensus assembly

optional arguments:
  -h, --help           show this help message and exit
  --min_cov INT        Minimum required coverage to not mask.
  --min_sgene_cov INT  The minimum required coverage of S gene (Spike) region.
  --width INT          Maximum line length for FASTA output (Default: 60, use 0 to print to a single line)
  --outdir STR         Directory to write output. (Default ./)
  --debug              Print full notes for every position (Default: only print changes
  --version            show program's version number and exit
```

#### `--min_cov` 
`--min_cov` sets the minimum per-base coverage to mask a sequence. For example, if `--min_cov` is set to 100 and a position has a coverage of only 80x, the base at that position will be set to *N*. This essentially allows you to take consensus assemblies from multiple workflows and uniformly apply coverage masking.

#### `--min_sgene_cov`
The *S gene* is important for lineage calls and therefore has its own parameter for coverage cutoffs. Please note though, the greater of `--min_cov` and `--min_sgene_cov` will be used for *S gene* coverage cutoffs.

#### `--width`
By default each line in the output FASTA file will be a maximum of 60 characters. `--width` will allow you to change this to your liking. Setting the value to `0` will put the sequence on a single line. Printing to a single line can be quite useful when you are planning to concatenate multiple consensus assemblies into a single file.

#### `--debug`
The `--debug` option will print out detail for all positions no matter if they were masked or not. This can be useful if you would like to manually review the coverages and basecalls across each position.


# Output Files
*rhea* outputs two files: a file containing details about any masking, and the final consensus assembly.


### `*-rhea.txt`
The `txt` file output by *rhea* includes information about any maskings that might have occured. It has 5 columns:

1. `position` - the position in the consensus assembly
2. `coverage` - the coverage observed at the given position
3. `observed_base` - the base observed in the original input consensus assembly
4. `masked_base` - the base that was included in the final masked consensus assembly
5. `note` - any details about why a position was masked


#### Example

```
position        coverage        observed_base   masked_base     note
21563   90     A       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21564   90     T       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21565   90     G       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21566   90     T       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21567   90     T       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21568   90     T       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
21569   90     G       N       Position did not meet the S gene coverage requirement, masking. (Observed 90x, Expected 100x)
```

You can use `--debug` to print out details for every position.

### `*-rhea.fasta`
The consensus assembly with positions masked by *rhea*.

```
><SAMPLE_NAME> SARS-CoV-2 consensus assembly with masking by Rhea [assembly_accession=MN908947.3] [length=29903]
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNAGATCT
GTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACT
CACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATC
TTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTT
TGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAAC
ACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGG
AGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGG
CTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAA
ACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACT
```


# Naming
Although *rhea* is not specific to a particular workflow, it was developed to supplement our *[Titan](https://github.com/theiagen/public_health_viral_genomics)* workflow. To keep the theme going, we chose to name this tool after the titan [Rhea](https://en.wikipedia.org/wiki/Rhea_(mythology)) which is also another one of [Saturn's Moons](https://en.wikipedia.org/wiki/Rhea_(moon)).
