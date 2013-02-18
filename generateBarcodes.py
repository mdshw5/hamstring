#!/usr/bin/env python

import os
import argparse
import hamstring

def main():
    parser = argparse.ArgumentParser(description='Generate Hamming DNA barcodes')
    parser.add_argument('out', help='output barcode file name')
    args = parser.parse_args()
    f = open(args.out, 'w')
    f.write('index base4 nucleotide gc\n')
    for i in range(0,256): ## range of decimal integers 4 ** number of data bits
        x = hamstring.base4Encode(i,4) ## generate base4 list 
        y = hamstring.generateHamming(x,3) ## generate hamming barcode with 3 parity bits
        z = ' '.join([str(y[a]) for a in y.keys()]) ## format results
        f.write('{0} '.format(str(i)))
        f.write(z)
        f.write('\n')
    f.close()

if __name__ == "__main__":
    main()
