#!/usr/bin/env python

import os
import argparse
import hamstring

def main():
    """ Read in a list of Hamming barcodes, one line at a time, perform checksum,
    and return corrected sequence, as well as error message to stdout
    """
    parser = argparse.ArgumentParser(description='Checksum a list of Hamming DNA barcodes')
    parser.add_argument('list', help='list of barcodes to check, one per line')
    args = parser.parse_args()
    with open(args.list,'r') as f:
        print 'in\tfixed\tchecksum'
        for line in f:
            a = line.rstrip()
            x = hamstring.decodeHamming(a,3)
            z = '\t'.join([a, x['nucleotide'], x['chksum']])
            print z
    f.close()

if __name__ == "__main__":
    main()
