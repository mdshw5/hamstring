# Readme

## Background

This python module defines functions that will generate quaternary Hamming barcodes, as well as perform checksum and error correction on supplied barcodes. The theory for generating quaternary encoded DNA Hamming barcodes comes from the publication [Bystrykh, L. V. (2012). Generalized DNA Barcode Design Based on Hamming Codes. PLoS ONE](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0036852).

## Usage

Generate Hamming DNA barcodes

    generateBarcodes.py [-h] out

    arguments:
      out         output barcode file name
      -h, --help  show this help message and exit

Checksum a list of Hamming DNA barcodes

    checkBarcodes.py [-h] list

    arguments:
      list        list of barcodes to check, one per line
      -h, --help  show this help message and exit

## Hamstring Module

`base4Encode(n,d)` is used to convert decimal notation *n* to quaternary notation with *d* leading digits. *example*: 

    base4Encode(22, 4)
    [0, 1, 1, 2]	    

`generateHamming(data,parity)` is used to generate DNA quaternary Hamming codes from list of quaternary digits *data* with *parity* number of parity bits.
*example*:

    generateHamming([0,1,1,2],3)
    {'parity': [1, 1, 0], 'nucleotide': 'CCAACCG', 'data': [0, 1, 1, 2], 'base4': '1100112'}

`decodeHamming(barcode,parity)` is used to decode *barcode* nucleotide Hamming string with *parity* number of parity bits, and perform error correction if needed.
*example*:

    decodeHamming('CCAACCG',3)
    {'nucleotide': 'CCAACCG', 'chksum': 'ok'}

    decodeHamming('CCATCCG',3)
    {'nucleotide': 'CCAACCG', 'chksum': 'T to A at 4'}

## Author Information

Matt Shirley - mdshw5'at'gmail'.'com - [http://mattshirley.com](http://mattshirley.com)

## To Do

- Implement binary to quaternary conversion (Done)
- Implement quaternary Hamming encoding function (Done)
- Implement quaternary Hamming decoding and checksumming function (Done)
- Implement quaternary Hamming distance calculation function
- Implement GC content calculations (Done)
- Implement barcode list checking (Done)
- Implement sequence redundancy calculations
- Implement fastq barcode error correction function