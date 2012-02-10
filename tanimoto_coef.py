#!/usr/bin/env python

from bitarray import bitarray

def counttc(refbit, tgtbit):
    """Measure Tanimoto coefficient of reference bit (refbit) and target bit (tgtbit)"""
    refcopy = refbit.copy()
    tgtcopy = tgtbit.copy()
    a = b = c = 0
    
    while refbit.length():
        if refbit.pop():
            a += 1
        if tgtbit.pop():
            b += 1
        if refcopy.pop() & tgtcopy.pop():
            c += 1
        if refcopy.length() < tgtcopy.length():
            tgtcopy.pop()

    tc = float(c)/(a+b-c)
    
    return tc

def gettcfromdict(dict1, dict2, residuechoice):
    bit1 = bitarray()
    bit2 = bitarray()
    for residue in residuechoice:
        bit1.append(dict1[residue])
        bit2.append(dict2[residue])
    tcfromdict = counttc(bit1, bit2)
    
    return tcfromdict

if __name__ == "__main__":
    print "This code is part of PyPLIF to measure Tanimoto coefficient from reference bit and target bit"
