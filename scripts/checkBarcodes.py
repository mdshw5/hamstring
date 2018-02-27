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
import hamstring


def main():
    """ Read in a list of Hamming barcodes, one line at a time, perform checksum,
    and return corrected sequence, as well as error message to stdout
    """
    parser = argparse.ArgumentParser(
        description='Checksum a list of Hamming DNA barcodes')
    parser.add_argument('list', help='list of barcodes to check, one per line')
    parser.add_argument(
        '-p',
        '--parity',
        type=int,
        default=3,
        help=
        'length of the parity bit e.g. 4 for Hamming8,4. default=%(default)s')
    args = parser.parse_args()
    with open(args.list, 'rU') as f:
        print('in\tfixed\tchecksum')  ## print header
        for line in f:
            a = line.rstrip()  ## remove newlines
            x = hamstring.decodeHamming(a, args.parity)
            z = '\t'.join([a, x['nucleotide'], x['chksum']])
            print(z)
    f.close()


if __name__ == "__main__":
    main()
