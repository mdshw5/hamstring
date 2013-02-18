#!/usr/bin/env python

import os
import argparse
import hamstring

def main():
    parser = argparse.ArgumentParser(description='Generate Hamming DNA barcodes')
    parser.add_argument('out', help='output barcode file name')
    args = parser.parse_args()
    f = open(args.out, 'w')
    for i in range(0,256):
        x = hamstring.base4Encode(i,4)
        y = hamstring.generateHamming(x,3)
        z = ' '.join([str(y[a]) for a in y.keys()])
        f.write('{0} '.format(str(i)))
        f.write(z)
        f.write('\n')
    f.close()

if __name__ == "__main__":
    main()
