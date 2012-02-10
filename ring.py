#!/usr/bin/env python

from numpy import cross

def getringdict(protein):
    rings = protein.GetSSSR()
    ringdict = {}
    ringnum  = 0
    for ring in rings:
        if ring.IsAromatic():
            ringdict[ringnum] = [protein.GetAtom(ring._path[0]).GetResidue().GetName(), ring._path]
            ringnum += 1
     
    return ringdict


def getringcross(coords):
    a = coords[0]-coords[1]
    b = coords[0]-coords[2]
    crossprod = cross(a,b)
    
    return crossprod


#  Suppose to be the distance between ring geometric center
#  but it doesn't work well, also the FingerprintLib by
#  Marcou & Rognan using the ring member coordinates instead
#  of the ring geometric center.
def ringdistance(ring1, ring2, mol1, mol2):
    for atom1 in ring1._path:
        for atom2 in ring2._path:
            atomdistance = mol1.GetAtom(atom1).GetDistance(mol2.GetAtom(atom2))
            if atomdistance <= 4.0:
                return True

    return False
