#! /usr/bin/env python3
"""
usage: rhea [-h] [--min_cov INT] [--min_spike_cov INT] [--width INT] [--version] SAMPLE FASTA VCF BAM

rhea (v0.0.1) - Snippy consensus (subs) with coverage masking.

positional arguments:
  SAMPLE               Sample name to use for naming outputs
  FASTA                Consensus assembly in FASTA format
  VCF                  VCF file from the consensus assembly
  BAM                  Per-base coverage of alignment

optional arguments:
  -h, --help           show this help message and exit
  --min_cov INT        Minimum required coverage to not mask.
  --min_spike_cov INT  For SARS-CoV-2 genomes, the minimum required coverage of Spike region.
  --width INT          Maximum line length for FASTA output (Default: 80, use 0 to print to a single line)
  --version            show program's version number and exit
"""
PROGRAM = "rhea"
VERSION = "0.0.1"
# Taken from https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3
SC2_GENOME_SIZE=29903
SPIKE_START=21563
SPIKE_STOP=25384
import sys
from collections import OrderedDict

def default_counts(name, length):
    counts = OrderedDict()
    for i in range(length):
        counts[str(i + 1)] = 0
    return {name: counts}

def get_per_base_coverage(coverage):
    """Extract per-base coverage from the BAM file"""
    import pysam
    bamfile = pysam.AlignmentFile(coverage, 'rb')
    reference = bamfile.references[0]
    coverages = default_counts(reference, bamfile.get_reference_length(reference))
    for row in bamfile.pileup(max_depth=0, stepper='nofilter'):
        pos = str(row.reference_pos + 1)
        depth = row.nsegments
        coverages[reference][pos] = depth
    return coverages

def read_vcf(vcf):
    """Get positions with a substitution."""
    subs = {}
    with open(vcf, 'rt') as vcf_fh:
        for line in vcf_fh:
            if not line.startswith("#"):
                line = line.split('\t')
                # 0 = accession, 1 = position
                if line[0] not in subs:
                    subs[line[0]] = {}
                subs[line[0]][line[1]] = True
    return subs


def read_fasta(fasta):
    """Parse the input FASTA file."""
    from Bio import SeqIO
    seqs = {}
    with open(fasta, 'r') as fasta_fh:
        for record in SeqIO.parse(fasta_fh,'fasta'):
            seqs[record.name] = str(record.seq)
    return seqs


def mask_sequence(sequence, coverages, subs, mincov):
    """Mask positions with low or no coverage in the input FASTA."""
    masked_seqs = {}
    
    for reference, vals in coverages.items():
        bases = []
        for pos, cov in vals.items():
            if cov >= mincov:
                # Passes
                if reference in subs:
                    if str(pos) in subs[reference]:
                        # Substitution
                        bases.append(sequence[reference][pos].lower())
                    else:
                        # Same as reference
                        bases.append(sequence[reference][pos])
                else:
                    # No SNPs, Same as reference
                    bases.append(sequence[reference][pos])
            elif cov:
                # Low coverage
                bases.append("N")
            else:
                # 0 coverage
                bases.append('n')

        if len(bases) != len(sequence[reference]):
            print(f'Masked sequence ({len(bases)} for {reference} not expected length ({len(sequence[reference])}).',
                file=sys.stderr)
            sys.exit(1)
        else:
            masked_seqs[reference] = bases

    return masked_seqs


def format_header(sample, reference, accession, length):
    """Return a newly formatted header."""
    title = f'Pseudo-seq with called substitutions and low coverage masked'
    return f'>gnl|{accession}|{sample} {title} [assembly_accession={reference}] [length={length}]'


def chunks(s, n):
    """
    Produce `n`-character chunks from `s`.
    https://stackoverflow.com/questions/7111068/split-string-by-count-of-characters
    """
    for start in range(0, len(s), n):
        yield s[start:start+n]


if __name__ == '__main__':
    import argparse as ap
    import sys

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - Snippy consensus (subs) with coverage masking.'
        )
    )
    parser.add_argument('sample', metavar="SAMPLE", type=str,
                        help='Sample name to use for naming outputs')
    parser.add_argument('fasta', metavar="FASTA", type=str,
                        help='Consensus assembly in FASTA format')
    parser.add_argument('vcf', metavar="VCF", type=str,
                        help='VCF file from the consensus assembly')
    parser.add_argument('bam', metavar="BAM", type=str,
                        help='Per-base coverage of alignment')
    parser.add_argument('--min_cov', metavar='INT', type=int, default=50,
                        help='Minimum required coverage to not mask.')
    parser.add_argument('--min_spike_cov', metavar='INT', type=int, default=100,
                        help='For SARS-CoV-2 genomes, the minimum required coverage of Spike region.')
    parser.add_argument('--width', metavar='INT', type=int, default=80,
                        help='Maximum line length for FASTA output (Default: 80, use 0 to print to a single line)')            
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    coverages = get_per_base_coverage(args.bam)
    sub_positions = read_vcf(args.vcf)
    seqs = read_fasta(args.fasta)
    masked_seqs = mask_sequence(seqs, coverages, sub_positions, args.min_cov)
    for accession, seq in masked_seqs.items():
        header = format_header(args.sample, args.reference, accession, len(seq))
        print(header)
        if args.width:
            for chunk in chunks(seq, args.width):
                print("".join(chunk))
        else:
            print(seq)
