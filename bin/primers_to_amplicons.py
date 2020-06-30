#!/usr/bin/env python
'''
Convert the nCoV primer scheme to a unique amplicon BED file.
'''

import sys
import argparse
import ncov.utils.primers as utils

parser = argparse.ArgumentParser(description='Create unique amplicon BED file')
parser.add_argument('-p', '--primers', help='Primer scheme in BED format')
parser.add_argument('--offset', default=30, help='Primer offset for coordinates')
parser.add_argument('-o', '--output', default='out.bed', help='filename to write BED to')
parser.add_argument('--pattern', default='nCoV-2019_', help='amplicon name pattern')
if len(sys.argv) <= 1:
    parser.print_help(sys.stderr)
    sys.exit(1)
args = parser.parse_args()

primers = utils.import_bed_file(bed=args.primers)
primer_pairs = utils.create_primer_pairs(primers=primers)
amplicons = utils.create_amplicon_range(primer_pairs=primer_pairs, pattern=args.pattern)
unique_amplicons = utils.create_unique_amplicons(amplicons, offset=args.offset)
with open(args.output, 'w') as file_o:
    for line in unique_amplicons:
        file_o.write('\t'.join([line]))
        file_o.write('\n')
file_o.close()
