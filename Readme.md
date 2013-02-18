# Readme

## Background

This python module defines functions that will generate quaternary Hamming barcodes, as well as perform checksum and error correction on supplied barcodes. The theory for generating quaternary encoded DNA Hamming barcodes comes from the publication [Bystrykh, L. V. (2012). Generalized DNA Barcode Design Based on Hamming Codes. PLoS ONE](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0036852). These two functions work very well. 

base4Encode(n,d) is used to convert decimal notation to quaternary notation.

generateHamming(data,parity) is used to generate DNA quaternary Hamming codes.

decodeHamming(barcode,parity) is used to decode nucleotide Hamming barcode sequence and perform error correction if needed.
 
