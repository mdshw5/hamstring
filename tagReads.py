#!/usr/bin/env python                                                                                       
import os
import argparse
import hamstring
import random
from Bio import SeqIO, Seq

def main():
    """ Read a fastq file and prepend a barcode to each read. Mutate the barcode e % of the time
    for nb barcodes.
    e = percent single error rate in barcode
    nb = number of barcodes in experiment
    TODO: implement barcode mutation and better fastq output
    """
    parser = argparse.ArgumentParser(description='Tag fastq reads with a barcode')
#    parser.add_argument('e', help='error rate for single barcode base errors. e.g. 0.10')
    parser.add_argument('nb', help='number of barcodes to generate')
    parser.add_argument('fastq', help='fastq file to process')
    parser.add_argument('out', help='name for new fastq file')
    args = parser.parse_args()
    ## generate barcodes 
    rNum = random.sample(range(0,256),int(args.nb))
    b4 = [hamstring.base4Encode(decimal,4) for decimal in rNum]
    barcodes = [hamstring.generateHamming(data,3)['nucleotide'] for data in b4]        
    ## open the fastq file and output file; read through, prepending barcodes
    o = open(args.out, 'w')
    with open(args.fastq, 'rU') as fq:
        for record in SeqIO.parse(fq, 'fastq'):
            newSeq = random.choice(barcodes) + record.seq
            qual = map(int,record.letter_annotations['phred_quality'])
            qual33 = [ascii + 33 for ascii in qual]
            qualText = ''.join([chr(qscore) for qscore in qual33])
            newQual = 'H'*7 + qualText
            ## really need to do something better than this for outfile
            o.write('@{0}'.format(record.id))
            o.write('\n')
            o.write(str(newSeq))
            o.write('\n')
            o.write('+')
            o.write('\n')
            o.write(str(newQual))
            o.write('\n')
    fq.close()
    o.close()
            

if __name__ == "__main__":
    main()
