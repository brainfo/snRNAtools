''' Map histone ChIP-seq reads '''

import argparse

parser = argparse.ArgumentParser(description=__doc__,epilog=EPILOG,
formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument(
    '--reference',
    help='Reference to map to'
)
parser.add_argument(
    '--chrom_sizes',
    help='chrom.sizes file for bedToBigBed'
)

args = parser.parse_args()

def align():

    cmd = 'bowtie2 -p 12 -q --local \
    -x '
