#!/usr/bin/env python

from bitarray import bitarray

def counttc(refbit, tgtbit):
	"""Measure Tanimoto coefficient of reference bit (refbit) and target bit (tgtbit)"""
	a = refbit.count()
	b = tgtbit.count()
	c = 0
	while refbit.length():
		if refbit.pop() & tgtbit.pop():
			c += 1
		if refbit.length() < tgtbit.length():
			tgtbit.pop()

	tc = float(c)/(a+b-c)
	
	return tc

def gettcfromdict(dict1, dict2, residuechoice):
	"""Return the Tanimoto coefficient from the given dictionary"""
	bit1 = bitarray()
	bit2 = bitarray()
	for residue in residuechoice:
		bit1.extend(dict1[residue])
		bit2.extend(dict2[residue])

	tcfromdict = counttc(bit1, bit2)
	
	return tcfromdict

def collectbit(dict1, residuechoice):
	"""This function will return bit array in string form"""
	bit = bitarray()
	for residue in residuechoice:
		bit.extend(dict1[residue])

	stringbit = str(bit).split("'")[1]

	return stringbit

if __name__ == "__main__":
	print "This code is part of PyPLIF to measure Tanimoto coefficient from reference bit and target bit"
