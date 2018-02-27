#!/usr/bin/env python
#
# Copyright (C) 2013 Matt Shirley
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import argparse
import random
from hamstring import *


def main():
    """ Read a fastq file and prepend a barcode to each read. Mutate the barcode e % of the time
    for nb barcodes.
    e = percent single error rate in barcode
    nb = number of barcodes in experiment
    """
    parser = argparse.ArgumentParser(
        description='Tag fastq reads with a barcode')
    parser.add_argument('nb', type=int, help='number of barcodes to generate')
    parser.add_argument('fastq', type=str, help='fastq file to process')
    parser.add_argument('out', type=str, help='name for new fastq file')
    parser.add_argument(
        '-e',
        '--erate',
        type=float,
        default=0.05,
        help='error rate for single barcode base errors. default=0.05')
    parser.add_argument(
        '-p',
        '--parity',
        type=int,
        default=3,
        help=
        'length of the parity bit e.g. 4 for Hamming8,4. default=%(default)s')
    args = parser.parse_args()
    ## generate barcodes
    rn = random.sample(range(0, 256), int(args.nb))
    b4 = [base4Encode(decimal, 4) for decimal in rn]
    barcodes = [
        generateHamming(data, args.parity)['nucleotide'] for data in b4
    ]
    print('Barcodes:\n' + '\n'.join(barcodes))
    ## generate weighted distribution
    wd = [True] * int(100 * args.e) + [False] * (100 - (100 * int(args.e)))
    ## open the fastq file and output file; read through, prepending barcodes
    with fastqReader(args.fastq) as fq, fastqWriter(args.out) as out:
        for record in fq:
            barcode = random.choice(barcodes)
            if random.choice(wd) is True:
                pos = random.choice([0, 1, 2, 3, 4, 5, 6])
                bases = ['A', 'C', 'G', 'T']
                bases = [b for b in bases if not b.endswith(barcode[pos])]
                barcode = list(barcode)
                barcode[pos] = random.choice(bases)
                barcode = ''.join(barcode)
            record.seq = barcode + record.seq
            record.qual = 'H' * 7 + record.qual
            out.write(record)


if __name__ == "__main__":
    main()
