#!/usr/bin/env python

import getopt, sys, os
from openbabel import OBMol, OBConversion
from optparse import OptionParser
from time import time
from glob import glob
from tanimoto_coef import *
from ring import *
from interactions import *


if __name__ == "__main__":
	x = time()
    

	#  Default configuration
	config    = "config.txt"

	#  Identifiers initialization
	protein_reference = ""
	ligand_reference = ""
	residue_of_choice = []
	protein_ligand_folder = ""
	output_file = ""

	parser = OptionParser()
	parser.add_option("-c", "--config", dest="config",
					  help="read config from FILE", metavar="FILE")
	parser.add_option("-o", "--output", dest="output_file",
					  help="write result to FILE", metavar="FILE")
	try:
		configread  = open(config, 'r')
	except:
		print "The config file: '", config, "' can not be found"
		sys.exit(1)

	configlines = [line for line in configread]
	configread.close()

	for line in configlines:
		options = line.split()
		if not options:
			continue
		if options[0] == "protein_reference":
			protein_reference = options[1]
		if options[0] == "ligand_reference":
			ligand_reference = options[1]
		if options[0] == "residue_of_choice":
			residue_of_choice = options[1:]
		if options[0] == "protein_ligand_folder":
			protein_ligand_folder = options[1]
		if options[0] == "output_file":
			output_file = options[1]

	try:
		os.chdir(protein_ligand_folder)
		mollisttemp = glob('*conf_01.mol2')
		mollist = [mol.split('_')[0] for mol in mollisttemp]
		mollist.sort()
		os.chdir('..')
		pbf = protein_ligand_folder + '/protein_bindingsite_fixed.mol2'
	except:
		print 'The protein ligand folder can not be found'
		sys.exit(1)
            

	# opening the molecule files
	conv = OBConversion()
	conv.SetInFormat("mol2")

	protfix  = OBMol()
	protref  = OBMol()
	ligref   = OBMol()
	docklig  = OBMol()
	dockprot = OBMol()

	conv.ReadFile(protfix, pbf)
	conv.ReadFile(protref, protein_reference)
	conv.ReadFile(ligref, ligand_reference)
    
	refresdict    = getresiduedict(protref)
	refringdict   = getringdict(protref)


	ringinteraction(refresdict, refringdict, residue_of_choice, protref, ligref)
	otherinteractions(refresdict, residue_of_choice, protref, ligref)


	conf_number = len(glob(protein_ligand_folder+'/'+mollist[0]+'*protein.mol2')) 

	cvsoutdict = {}
	bitarraydict = {}
	for compound in mollist:
		cvsoutdict[compound] = []
		bitarraydict[compound] = []
		for conf in range(conf_number):
			conf += 1
			fixringdict   = getringdict(protfix)
			fixresdict    = getresiduedict(protfix)

			base_name = protein_ligand_folder + '/' + compound + '_entry_00001_conf_' + str(conf).zfill(2)
			ligand_file = base_name + '.mol2'
			protein_file = base_name + '_protein.mol2'
			conv.ReadFile(docklig, ligand_file)
			conv.ReadFile(dockprot, protein_file)
			ringinteraction(fixresdict, fixringdict, residue_of_choice, protfix, docklig)
			otherinteractions(fixresdict, residue_of_choice, protfix, docklig)
			hbonddockprot(fixresdict, residue_of_choice, dockprot, docklig)

			tc = gettcfromdict(refresdict, fixresdict, residue_of_choice)
			stringbit = collectbit(fixresdict, residue_of_choice)
			cvsoutdict[compound].append(tc)
			bitarraydict[compound].append(stringbit)
    
	refstringbit = collectbit(refresdict, residue_of_choice)
	
	outfile = open(output_file, 'w')
	outfile.write("                   ")
	for res in residue_of_choice:
		outfile.write(res.ljust(7))
	outfile.write('\n')
	outfile.write(ligref.GetTitle().ljust(18))
	outfile.write(" %s" % refstringbit)
	outfile.write('\n')
	for compound in cvsoutdict:
		for stringbit, tc in zip(bitarraydict[compound], cvsoutdict[compound]):
			outfile.write(compound.ljust(18))
			outfile.write(" %s" % stringbit)
			outfile.write(" %.3f" % tc)
			outfile.write('\n')

#        for tc in cvsoutdict[compound]:
                    
	outfile.close()
	y = time()
	print 'Total time taken %.3f s.' % (y-x)





