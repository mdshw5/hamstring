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
from hamstring import *


def main():
    """ Read fastq files from input file, parse barcode, calculate checksum, and write
    a new, corrected fastq file as well as a log of operations. Note that the barcode is
    expected to be the first n nucleotides from the 5' end.

    """
    parser = argparse.ArgumentParser(
        description='Check and fix barcodes in fastq file')
    parser.add_argument(
        'list', help='list of barcodes used in experiment, one per line')
    parser.add_argument('fastq', help='fastq file to process')
    parser.add_argument('out', help='name for new fastq file')
    parser.add_argument(
        '-s',
        '--strict',
        action="store_true",
        help='change all barcodes not in list to \'N\'')
    parser.add_argument(
        '-p',
        '--parity',
        type=int,
        default=3,
        help=
        'length of the parity bit e.g. 4 for Hamming8,4. default=%(default)s')
    args = parser.parse_args()
    n = 7  ## nucleotide length of barcode sequence
    ## open the list of barcodes and read them into a list object
    with open(args.list, 'rU') as f:
        barcodes = []
        for line in f:
            barcodes.append(line.rstrip())
    f.close()
    ## open the fastq file and output file; read through, performing checksums and write corrected reads
    with fastqReader(args.fastq) as fq, fastqWriter(args.out) as out:
        for record in fq:
            barcode = record.seq[:n]
            decode = decodeHamming(barcode, args.parity)
            seq = list(record.seq)
            if (decode['chksum'] != 'ok'
                    and any(decode['nucleotide'] in s for s in barcodes)):
                seq[:n] = list(decode['nucleotide'])
                messg = 'corrected ' + decode['chksum'] + ' in read ' + record.name
                print(messg)
            elif (args.strict and decode['chksum'] != 'ok'
                  and not any(decode['nucleotide'] in s for s in barcodes)):
                seq[:n] = list('N' * len(barcode))
                messg = 'discarded barcode ' + barcode + ' in read ' + record.name
                print(messg)
            record.seq = ''.join(seq)
            out.write(record)


if __name__ == "__main__":
    main()
