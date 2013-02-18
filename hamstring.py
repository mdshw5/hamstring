#!/usr/bin/env python
# Principles used in this script are adapted from the 2012 publication:
# Bystrykh, L. V. (2012). Generalized DNA Barcode Design Based on Hamming Codes. 
# PLoS ONE, 7(5), e36852. doi:10.1371/journal.pone.0036852.t004

from __future__ import division

def base4Encode(n,d):
    """Convert decimal notation to quaternary notation
    We will use division and modulus recursively

    n = decimal number
    d = number of digits for quaternary representation

    """
    alphabet = [0,1,2,3]
    quat = []
    base = len(alphabet)
    while d: ## recursively calculate modulus
        remainder = n % base
        n = int(n / base)
        quat.append(alphabet[remainder])
        d-=1
    quat.reverse()
    return quat
  
def generateHamming(data,parity):
    """Generate quaternary Hamming codes

    data = quaternary number list
    parity = number of parity bits to implement

    """
    d = len(data) ## number of coding bits
    N = ['A','C','G','T'] ## nucleotide dictionary
    # Generate barcodes
    p1 = (d - sum([data[i] for i in [0,1,3]]) % d) % d
    p2 = (d - sum([data[i] for i in [0,2,3]]) % d) % d
    p3 = (d - sum([data[i] for i in [1,2,3]]) % d) % d
    h4 = [p1,p2,data[0],p3,data[1],data[2],data[3]]
    ## substitute nucleotide for quaternary encoding
    hN = map(lambda x: N[x], h4)
    s4 = ''.join([str(i) for i in h4])
    sN = ''.join(hN)
    gc = percentGC(sN)
    z = {'gc':gc, 'base4':s4, 'nucleotide':sN}
    return z

def percentGC(x):
    """ Calculate the percent GC content of a nucleotide string
    Return result rounded to two decimal places

    x = nucleotide string
    """
    y = list(x)
    l = len(y)
    Q = {'A':0, 'C':1, 'G':1, 'T':0}
    z = sum(map(lambda x: Q[x], y)) / l
    return round(z,2)

def smashBase(x):
    """ Smash a base4 number to base2 
    
    x = a base4 encoded value
    """
    if x > 0:
        x = 1
    return x

def decodeHamming(barcode,parity):
    """ Decode nucleotide Hamming barcode sequence and perform error correction

    """
    d = len(barcode) - parity
    hN = list(barcode)
    Q = {'A':0, 'C':1, 'G':2, 'T':3}
    N = {0:'A', 1:'C', 2:'G', 3:'T'}
    h4 = map(lambda x: Q[x], hN)
    p1 = sum([h4[i] for i in [0,2,4,6]]) % d
    p2 = sum([h4[i] for i in [1,2,5,6]]) % d
    p3 = sum([h4[i] for i in [3,4,5,6]]) % d
    errType = max(p1,p2,p3) ## determine the type of error for later correction
    ## determine position of error
    binErr = map(smashBase, [p3,p2,p1]) ## get the reversed "binary" version of parity bits
    errPos = int(''.join([str(i) for i in binErr]),2) ## error position
    if errType != 0:
        sFalse = h4[errPos -1]
        sTrue = h4[errPos - 1] - errType
        if sTrue < 0:
            sTrue+=4
        h4[errPos - 1] = sTrue
        hN = map(lambda x: N[x], h4)
        sN = ''.join([str(i) for i in hN])
        errMsg = ' '.join([N[sTrue],'>',N[sFalse],'at pos',str(errPos)])
        return {'nucleotide':sN, 'chksum':errMsg}
    else: 
        return {'nucleotide':barcode, 'chksum':'ok'}
