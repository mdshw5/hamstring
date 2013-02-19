#!/usr/bin/env python                                                                                       
import os
import argparse
import hamstring
from Bio import SeqIO

def main():
    """ Read fastq files from input file, parse barcode, calculate checksum, and write 
    a new, corrected fastq file as well as a log of operations. Note that the barcode is 
    expected to be the first n nucleotides from the 5' end.

    """
    parser = argparse.ArgumentParser(description='Check and fix barcodes in fastq file')
    parser.add_argument('list', help='list of barcodes used in experiment, one per line')
    parser.add_argument('fastq', help='fastq file to process')
    parser.add_argument('out', help='name for new fastq file')
    args = parser.parse_args()
    n = 7 ## nucleotide length of barcode sequence
    ## open the list of barcodes and read them into a list object
    with open(args.list,'rU') as f:
        barcodes = []
        for line in f:
            barcodes.append(line.rstrip())
    f.close()
    ## open the fastq file and output file; read through, performing checksums and write corrected reads
    o = open(args.out, 'w')
    with open(args.fastq, 'rU') as fq:
        for record in SeqIO.parse(fq, 'fastq'):
            seq = list(record.seq)
            qual = map(int,record.letter_annotations['phred_quality'])
            qual33 = [ascii + 33 for ascii in qual]
            qualPhred = ''.join([chr(qscore) for qscore in qual33])
            barcode = ''.join(seq[:n])
            decode = hamstring.decodeHamming(barcode,3)
            if (decode['chksum'] != 'ok' and any(decode['nucleotide'] in s for s in barcodes)):
                seq[:n] = list(decode['nucleotide'])
                messg = 'corrected ' + decode['chksum'] + ' in read ' + record.id
                print messg
            o.write('@{0}'.format(record.id))
            o.write('\n')
            o.write(''.join(seq))
            o.write('\n')
            o.write('+')
            o.write('\n')
            o.write(str(qualPhred))
            o.write('\n')
    fq.close()
    o.close()
            

if __name__ == "__main__":
    main()
