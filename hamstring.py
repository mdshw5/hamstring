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
#
# Principles used in this script are adapted from the 2012 publication:
# Bystrykh, L. V. (2012). Generalized DNA Barcode Design Based on Hamming Codes.
# PLoS ONE, 7(5), e36852. doi:10.1371/journal.pone.0036852.t004
#
#
# To Do
#- Implement binary to quaternary conversion (Done)
#- Implement quaternary Hamming encoding function (Done)
#- Implement quaternary Hamming decoding and checksumming function (Done)
#- Implement quaternary Hamming distance calculation function
#- Implement GC content calculations (Done)
#- Implement barcode list checking (Done)
#- Implement sequence redundancy calculations
#- Implement fastq barcode error correction function (Done)
#- Implement barcode mutation rate in simulation tool

from __future__ import division
from collections import namedtuple


class read(object):
    """
    A class to hold features from fastq reads.
    """

    def __init__(self, name, seq, strand, qual):
        self.name = name
        self.seq = seq
        self.strand = strand
        self.qual = qual
        self.dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

    def __getitem__(self, key):
        return self.__class__(self.name, self.seq[key], self.strand,
                              self.qual[key])

    def index(self):
        return list(range(len(self.seq)))

    def seqlen(self):
        return len(self.seq)

    def reverse(self):
        """ Reverse the order of self """
        return self.__class__(self.name, self.seq[::-1], self.strand,
                              self.qual[::-1])

    def complement(self):
        """ Take the compliment of read.seq """
        compseq = ''.join([self.dict[x] for x in self.seq])
        return self.__class__(self.name, compseq, self.strand, self.qual)

    def revcomplement(self):
        """ Take the reverse compliment of read.seq """
        revcompseq = ''.join([self.dict[x] for x in self.seq])[::-1]
        return self.__class__(self.name, revcompseq, self.strand,
                              self.qual[::-1])

    def trim3(self, start, end):
        """ Trim all read class elements from the 3' end to the start of the adapter sequence alignment """
        self.seq, self.qual = self.seq[:start], self.qual[:start]
        return self

    def trim5(self, start, end):
        """ Trim all read class elements from the 5' start to the end  of the adapter sequence alignment """
        self.seq, self.qual = self.seq[end:], self.qual[end:]
        return self

    def trim53(self, start, end):
        """ Trim all read class elements to the 1-based 'trim' values"""
        self.seq, self.qual = self.seq[start:end], self.qual[start:end]
        return self


class fastqReader:
    """
    A class to read the name, sequence, strand and qualities from a fastq file

    file = the file name of a fastq file
    """

    def __init__(self, file):
        self.file = open(file, 'rU')
        self.read = read

    def __iter__(self):
        """
        Return read class: (name, sequence, strand, qualities).
        """
        for i, line in enumerate(self.file):
            if i % 4 == 0:
                name = line.strip()[1:]
            elif i % 4 == 1:
                sequence = line.strip()
            elif i % 4 == 2:
                strand = line.strip()
            elif i % 4 == 3:
                qualities = line.rstrip('\n\r')
                yield self.read(name, sequence, strand, qualities)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class fastqWriter:
    """ Take a read class object and file name, open file and write read """

    def __init__(self, file):
        self.file = open(file, 'w')

    def write(self, read):
        self.file.write('@' + read.name.split('/')[0] + '\n')
        self.file.write(read.seq + '\n')
        self.file.write(read.strand + '\n')
        self.file.write(read.qual + '\n')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


def base4Encode(n, d):
    """ Convert decimal notation to quaternary notation
    We will use division and modulus recursively

    n = decimal number
    d = number of digits for quaternary representation

    >>> base4Encode(22, 4)
    [0, 1, 1, 2]
    """
    alphabet = [0, 1, 2, 3]
    quat = []
    base = len(alphabet)
    while d:  ## recursively calculate modulus
        remainder = n % base
        n = int(n / base)
        quat.append(alphabet[remainder])
        d -= 1
    quat.reverse()
    return quat


def generateHamming(data, parity):
    """ Generate quaternary Hamming codes
    data = quaternary number list
    parity = number of parity bits to implement

    >>> generateHamming([0,1,1,2], 3)
    Barcode(base4='1100112', nucleotide='CCAACCG', gc=0.71)

    >>> generateHamming([0,1,1,2], 4)
    Barcode(base4='11001122', nucleotide='CCAACCGG', gc=0.75)
    """
    d = len(data)  ## number of coding bits
    if d != 4:
        raise ValueError("The number of data bits in the barcode exceeds 4.")
    N = ['A', 'C', 'G', 'T']  ## nucleotide dictionary
    # Generate barcodes
    p1 = (d - sum([data[i] for i in [0, 1, 3]]) % d) % d
    p2 = (d - sum([data[i] for i in [0, 2, 3]]) % d) % d
    p3 = (d - sum([data[i] for i in [1, 2, 3]]) % d) % d
    h4 = [p1, p2, data[0], p3, data[1], data[2], data[3]]
    if parity == 4:
        p4 = (d - sum(h4) % d) % d
        h4.append(p4)  # add extended parity bit
    elif parity != 3:
        raise ValueError("The parity argument must be 3 or 4.")
    ## substitute nucleotide for quaternary encoding
    hN = [N[x] for x in h4]
    s4 = ''.join([str(i) for i in h4])
    sN = ''.join(hN)
    gc = percentGC(sN)
    z = {'gc': gc, 'base4': s4, 'nucleotide': sN}
    return namedtuple('Barcode', z.keys())(**z)


def percentGC(x):
    """ Calculate the percent GC content of a nucleotide string
    Return result rounded to two decimal places

    x = nucleotide string
    """
    y = list(x)
    l = len(y)
    Q = {'A': 0, 'C': 1, 'G': 1, 'T': 0}
    z = sum([Q[x] for x in y]) / l
    return round(z, 2)


def smashBase(x):
    """ Smash a base4 number to base2

    x = a base4 encoded value
    """
    if x > 0:
        x = 1
    return x


def decodeHamming(barcode, parity):
    """ Decode nucleotide Hamming barcode sequence and perform error correction

    >>> decodeHamming('CCAACCG', 3)
    CheckedBarcode(nucleotide='CCAACCG', chksum='ok')

    >>> decodeHamming('CCAACCGG', 4)
    CheckedBarcode(nucleotide='CCAACCGG', chksum='ok')

    >>> decodeHamming('CCATCCG', 3)
    CheckedBarcode(nucleotide='CCAACCG', chksum='A > T at pos 4')

    >>> decodeHamming('CCATCCGG', 4)
    CheckedBarcode(nucleotide='CCAACCGG', chksum='A > T at pos 4')

    >>> decodeHamming('TCATCCGG', 4)
    CheckedBarcode(nucleotide='NNNNNNNN', chksum='bad')
    """
    d = len(barcode) - parity
    if d != 4:
        raise ValueError("The number of data bits in the barcode exceeds 4.")
    hN = list(barcode)
    Q = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    N = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    if any('N' in s for s in hN):
        z = {'nucleotide': str('N' * len(barcode)), 'chksum': 'bad'}
    else:
        h4 = list([Q[x] for x in hN])
        p1 = sum([h4[i] for i in [0, 2, 4, 6]]) % d
        p2 = sum([h4[i] for i in [1, 2, 5, 6]]) % d
        p3 = sum([h4[i] for i in [3, 4, 5, 6]]) % d
        errType = max(p1, p2,
                      p3)  ## determine the type of error for later correction
        ## determine position of error
        binErr = list(map(
            smashBase,
            [p3, p2, p1]))  ## get the reversed "binary" version of parity bits
        errPos = int(''.join([str(i) for i in binErr]), 2)  ## error position
        chksum = 'ok'
        if errType != 0:
            sFalse = h4[errPos - 1]
            sTrue = h4[errPos - 1] - errType
            if sTrue < 0:
                sTrue += 4
            h4[errPos - 1] = sTrue
            pp1 = sum([h4[i] for i in [0, 2, 4, 6]]) % d
            pp2 = sum([h4[i] for i in [1, 2, 5, 6]]) % d
            pp3 = sum([h4[i] for i in [3, 4, 5, 6]]) % d
            if parity == 4:
                pp4 = sum(h4) % d  # p4 checks the original Hamming7,4 code
            elif parity == 3:
                pp4 = 0
            elif parity != 3:
                raise ValueError("The parity argument must be 3 or 4.")
            if max(pp1, pp2, pp3, pp4) > 0:
                z = {'nucleotide': str('N' * len(barcode)), 'chksum': 'bad'}
            else:
                hN = [N[x] for x in h4]
                sN = ''.join([str(i) for i in hN])
                chksum = ' '.join(
                    [N[sTrue], '>', N[sFalse], 'at pos',
                     str(errPos)])
                z = {'nucleotide': sN, 'chksum': chksum}

        else:
            z = {'nucleotide': barcode, 'chksum': chksum}
    return namedtuple('CheckedBarcode', z.keys())(**z)


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
