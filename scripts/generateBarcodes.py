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

import argparse
import hamstring


def main():
    parser = argparse.ArgumentParser(
        description='Generate Hamming DNA barcodes')
    parser.add_argument('out', help='output barcode file name')
    parser.add_argument(
        '-p',
        '--parity',
        type=int,
        default=3,
        help=
        'length of the parity bit e.g. 4 for Hamming8,4. default=%(default)s')
    args = parser.parse_args()
    f = open(args.out, 'w')
    f.write('index base4 nucleotide gc\n')
    for i in range(0,
                   256):  ## range of decimal integers 4 ** number of data bits
        x = hamstring.base4Encode(i, 4)  ## generate base4 list
        y = hamstring.generateHamming(x, args.parity)
        z = ' '.join([str(y[a]) for a in list(y.keys())])  ## format results
        f.write('{0} '.format(str(i)))
        f.write(z)
        f.write('\n')
    f.close()


if __name__ == "__main__":
    main()
