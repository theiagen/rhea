#! /usr/bin/env python3
"""
usage: rhea [-h] [--min_cov INT] [--min_sgene_cov INT] [--width INT] [--outdir STR] [--debug] [--version] SAMPLE FASTA VCF BAM

rhea (v0.0.1) - SARS-CoV-2 consensus assembly quality checks and corrections

positional arguments:
  SAMPLE               Sample name to use for naming outputs
  FASTA                Consensus assembly in FASTA format
  VCF                  VCF file from the consensus assembly
  BAM                  Per-base coverage of alignment

optional arguments:
  -h, --help           show this help message and exit
  --min_cov INT        Minimum required coverage to not mask.
  --min_sgene_cov INT  The minimum required coverage of S gene (Spike) region.
  --width INT          Maximum line length for FASTA output (Default: 60, use 0 to print to a single line)
  --outdir STR         Directory to write output. (Default ./)
  --debug              Print full notes for every position (Default: only print changes
  --version            show program's version number and exit
"""
PROGRAM = "rhea"
VERSION = "0.0.1"

# Taken from https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3
SC2_REFERENCE="MN908947.3"
SC2_GENOME_SIZE=29903
SC2_SGENE_START=21563
SC2_SGENE_STOP=25384
import sys
from collections import OrderedDict

def default_counts(length: int) -> dict:
    """
    Generate a dictionary with each position set to 0 coverage

    Args:
        length (int): The number of positions to set to 0 coverage

    Returns:
        dict: Keys are reference position with the value set to 0 (e.g. dict['42'] -> 0)
    """
    counts = OrderedDict()
    for i in range(length):
        # range is 0-based (e.g. range(3) -> [0,1,2,3]), so add 1 to get reference position
        counts[str(i + 1)] = 0
    return counts

def read_bam(bam_file: str) -> dict:
    """
    Extract per-base coverage from the BAM file.

    Args:
        bam_file (string): Path to the input BAM file

    Returns:
        dict: Keys are reference position with the value set to observed coverage (e.g. dict['42'] -> 777)
    """
    import pysam
    bam = pysam.AlignmentFile(bam_file, 'rb')

    # `pileup()` only outputs positions with coverage, so preset all positions to 0 coverage
    coverages = default_counts(bam.get_reference_length(SC2_REFERENCE))

    # Docs: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignmentFile.pileup
    for pileupcolumn in bam.pileup(max_depth=0, stepper='nofilter'):
        # Update with observed coverage, 'reference_pos' is 0-based, so add 1
        # Docs: https://pysam.readthedocs.io/en/latest/api.html#pysam.PileupColumn
        coverages[str(pileupcolumn.reference_pos + 1)] = pileupcolumn.nsegments
    return coverages

def read_vcf(vcf_file: str) -> dict:
    """
    Read VCF file with observed mutations against the reference

    Args:
        vcf_file (string): Path to the input BAM file

    Returns:
        dict: Keys are reference position with the value set to VCF entry
    """
    import vcf
    records = OrderedDict()
    with open(vcf_file, 'r') as vcf_fh:
        for record in vcf.Reader(vcf_fh):
            records[str(record.POS)] = record
    return records

def read_fasta(fasta_file: str) -> dict:
    """
    Read FASTA formatted consensus assembly

    Args:
        fasta_file (string): Path to the input FASTA file

    Returns:
        dict: Keys are reference position with the value set to the observed nucleotide
    """
    from Bio import SeqIO
    seqs = []

    with open(fasta_file, 'r') as fasta_fh:
        for record in SeqIO.parse(fasta_fh,'fasta'):
            seqs.append(str(record.seq))

    # Run checks
    if len(seqs) == 0:
        # Expect single contig consensus assembly, found empty FASTA
        print(f"ERROR: Please verify input consensus assembly ({fasta_file}) is not empty", file=sys.stderr)
        sys.exit(2)
    elif len(seqs) > 1:
        # Expect single contig consensus assembly, found multiple FASTA entries
        print(f"ERROR: Input FASTA consensus assembly has multiple entries (found {len(seqs)} entries), expected only 1)",
              file=sys.stderr)
        sys.exit(2)
    elif len(seqs[0]) != SC2_GENOME_SIZE:
        # consensus assembly not equal to the expected length
        print(f"ERROR: Input consensus assembly has a length of {len(seqs[0])} bp, expected {SC2_GENOME_SIZE} bp",
              file=sys.stderr)
        sys.exit(2)
    else:
        # Basic checks passed, convert to per-base dict
        consensus = OrderedDict()
        for i, base in enumerate(seqs[0]):
            # enumerate is 0 based, add 1
            consensus[str(i + 1)] = base
        return consensus


def check_consensus(consensus: dict, coverages: dict, mutations: dict, min_cov: int, min_sgene_cov: int) -> dict:
    """
    Mask positions with low or no coverage in the input FASTA.

    Args:
        consensus (dict): Consensus assembly sequence
        coverages (dict): Per-base coverage of the consensus assembly
        mutations (dict): Mutations present in the input VCF file
        min_cov (int): The minimum coverage, coverages below will be masked.
        min_sgene_cove: The minimum coverage of S gene (Spike) region, coverages below will be masked.

    Returns:
        dict: Keys are reference positions with the value set to masked bases with notes
    """
    masked_consensus = OrderedDict()
    masked_sequence = []
    has_reference_inclusions = False
    for pos, base in consensus.items():
        cov = coverages[pos]
        passed = True

        if cov == 0 and base.lower() != 'n':
            # Reference base was included, warn the user
            passed = False
            masked_consensus[pos] = {
                'observed_base': base,
                'masked_base': 'N',
                'coverage': cov,
                'note': f'WARNING: Potential reference base inclusion. Observed {base} with 0x coverage at this position.',
                'passed': passed
            }
            print(f'WARNING: Potential reference base inclusion. Observed {base} with 0 coverage at position {pos}.', file=sys.stderr)
        elif cov > 0 and base.lower() == 'n':
            passed = False
            masked_consensus[pos] = {
                'observed_base': base,
                'masked_base': base,
                'coverage': cov,
                'note': f'WARNING: Potential reference base inclusion. Observed {base} with {cov}x coverage at this position.',
                'passed': passed
            }
        else:
            if int(pos) >= SC2_SGENE_START and int(pos) <= SC2_SGENE_STOP:
                # In the S gene region, use S gene-specific coverage cutoff, exclude bases already masked (N)
                if cov < min_sgene_cov and base.lower() != 'n':
                    passed = False
                    masked_consensus[pos] = {
                        'observed_base': base,
                        'masked_base': 'N',
                        'coverage': cov, 
                        'note': f'Position did not meet the S gene coverage requirement, masking. (Observed {cov}x, Expected {min_sgene_cov}x)',
                        'passed': passed
                    }
            else:
                # Use standard coverage cutoff, exclude bases already masked (N)
                if cov < min_cov and base.lower() != 'n':
                    passed = False
                    masked_consensus[pos] = {
                        'observed_base': base,
                        'masked_base': 'N',
                        'coverage': cov,
                        'note': f'Position did not meet the coverage requirement, masking. (Observed {cov}x, Expected {min_cov}x)',
                        'passed': passed
                    }
        
        if passed:
            # Everything passed, used the observed base
            masked_consensus[pos] = {
                'observed_base': base,
                'masked_base': base,
                'coverage': cov,
                'note': f'',
                'passed': passed
            }
        
        # Rebuild the consensus
        masked_sequence.append(masked_consensus[pos]['masked_base'])

    # Make sure consensus lengths match
    if len(masked_sequence) != SC2_GENOME_SIZE:
        # Masked consensus assembly not equal to the expected length
        print(f"ERROR: Masked consensus assembly has a length of {len(masked_sequence)} bp, expected {SC2_GENOME_SIZE} bp",
            file=sys.stderr)
        sys.exit(3)
    return [masked_consensus, "".join(masked_sequence)]


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
            f'{PROGRAM} (v{VERSION}) - SARS-CoV-2 consensus assembly quality checks and corrections'
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
    parser.add_argument('--min_sgene_cov', metavar='INT', type=int, default=100,
                        help='The minimum required coverage of S gene (Spike) region.')
    parser.add_argument('--width', metavar='INT', type=int, default=60,
                        help='Maximum line length for FASTA output (Default: 60, use 0 to print to a single line)')
    parser.add_argument('--outdir', metavar="STR", type=str, default="./",
                        help='Directory to write output. (Default ./)')
    parser.add_argument('--debug', action='store_true',
                        help='Print full notes for every position (Default: only print changes')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Read inputs: BAM file (per-base coverages), VCF (mutations), FASTA (consensus assembly)
    coverages = read_bam(args.bam)
    mutations = read_vcf(args.vcf)
    consensus = read_fasta(args.fasta)

    # Check the consensus, mask low coverages, warn reference inclusions
    masked_details, masked_consensus = check_consensus(consensus, coverages, mutations, args.min_cov, max(args.min_cov, args.min_sgene_cov))

    # Write any details regarding masking
    with open(f'{args.outdir}/{args.sample}-rhea.txt', 'wt') as rhea_txt:
        rhea_txt.write(f'position\tcoverage\tobserved_base\tmasked_base\tnote\n')
        for pos, vals in masked_details.items():
            if args.debug or vals['passed'] == False:
                rhea_txt.write(f'{pos}\t{vals["coverage"]}\t{vals["observed_base"]}\t{vals["masked_base"]}\t{vals["note"]}\n')

    # Write the final masked sequenced
    with open(f'{args.outdir}/{args.sample}-rhea.fasta', 'wt') as rhea_fasta:
        rhea_fasta.write(f'>{args.sample} SARS-CoV-2 consensus assembly with masking by Rhea [assembly_accession={SC2_REFERENCE}] [length={len(masked_consensus)}]\n')
        if args.width:
            for chunk in chunks(masked_consensus, args.width):
                seq = "".join(chunk)
                rhea_fasta.write(f'{seq}\n')
        else:
            rhea_fasta.write(f'{masked_consensus}\n')
