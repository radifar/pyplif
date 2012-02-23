#!/usr/bin/env python

from openbabel import OBResidueAtomIter


def setformalcharge(atom):
	#  This valence model is derived from OEChem valence model      #
	#  with slight modification. Atoms are assigned with formal     #
	#  charge based on their bond order. Nitrogen with bond order   #
	#  1 is given charge +1 since the output from PLANTS            #
	#  'protein_bindingsite_fixed.mol2' shows hydrogens detached    #
	#  from ammonium moeiety in Lysine. Whereas oxygen with bond    #
	#  order 1 but attached to hydrogen like the one appears in     #
	#  'docked_proteins.mol2' assigned with formal charge zero.     #
	#  																#
	#  OEChem Valence Model:										#
	#  http://www.eyesopen.com/docs/toolkits/current/html/OEChem_TK-c++/valence.html

	valence_model = {'boron':[0,2,1,0,-1,-2], 'nitrogen':[0,1,-1,0,1,0],
		            'oksigen':[0,-1,0,1,2,1], 'phosphorus':[0,4,3,2,1,0,-1,-2],
		            'sulfur':[0,-1,0,1,2,1,0,-1], 'chlor':[0,0,1,0,3]}

	#  Using the atomic number as the SMARTS pattern like so:       #
	#  #{atomic number}, so #5 = boron, #7 = nitrogen, and so on    #
	if atom.MatchesSMARTS('[#5,#7,#8,#15,#16,Cl]'):
		if atom.IsOxygen():
			atom.SetFormalCharge(valence_model['oksigen'][atom.BOSum()])
			if atom.BOSum() == 1:
				parentresidue = atom.GetResidue()
				for neighbor in OBResidueAtomIter(parentresidue):
					if atom.IsConnected(neighbor) & neighbor.IsHydrogen():
						atom.SetFormalCharge(0)
		elif atom.IsNitrogen():
			atom.SetFormalCharge(valence_model['nitrogen'][atom.BOSum()])
		elif atom.MatchesSMARTS('[#5]'):
			atom.SetFormalCharge(valence_model['boron'][atom.BOSum()])
		elif atom.IsPhosphorus():
			atom.SetFormalCharge(valence_model['phosphorus'][atom.BOSum()])
		elif atom.IsSulfur():
			atom.SetFormalCharge(valence_model['sulfur'][atom.BOSum()])
		elif atom.MatchesSMARTS('[Cl]'):
			atom.SetFormalCharge(valence_model['chlor'][atom.BOSum()])

def isnonpolar(atom):
	if atom.MatchesSMARTS('[#6,#16,F,Cl,Br,I]'):
		return True

	return False


if __name__ == "__main__":
	print "This code is part of PyPLIF to set atom flags based on their type."
