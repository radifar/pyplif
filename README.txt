PyPLIF is a program/script written in Python to analyze protein-ligand
interaction from the molecular docking result. It relies on OpenBabel
module to read and process the molecule files from the molecular
docking output.

The current version only support the Sybil mol2 format, but it is planned
to support other format in the future release such as pdbqt & dlg from
Autodock and Autodock Vina.

Unfortunately the current version is specially designed for analyzing the
docking results from PLANTS <http://www.tcd.uni-konstanz.de/research/plants.php>


Installation
------------

The instruction here is for installing PyPLIF as a module, for installing
PyPLIF as the stand-alone program check INSTALL.txt

To install PyPLIF module just use this command:
	sudo python setup.py install


Getting Started
---------------

Requirements:
* To use this program you should use the output from PLANTS.
* The ligand's conformation should be separated, that's the write_multi_mol2
should be disabled by putting this line into the PLANTS configuration file:
write_multi_mol2 0

Here's how to use PyPLIF step by step:

* Makes sure PyPLIF already installed, check INSTALL.txt for the complete
instruction
* Prepare the PyPLIF input by docking any ligand and protein you want
* Prepare config file for PyPLIF, the name can be anything but by default
PyPLIF will looking for 'config.txt', the config file example can be found in
docs folder. The content of it should be self-explaining.
* Also prepare the reference ligand and the reference protein conformation
in mol2 format.
* Run pyplif simply by entering the command 'pyplif'.
* After pyplif finished the calculation the output file will appear showing
the molecule file name, score, interaction bit, and Tc-IFP (Tanimoto coefficient
Interaction Fingerprinting)

Tips:
It is highly suggested to use the protein binding site instead of the whole
protein since using the binding site is much faster than using the whole protein.
You can use PLANTS to produce the binding site using the following command:
PLANTS --mode bind molecule.mol2 x protein.mol2
where the 'molecule.mol2' is the ligand file, x is the additional distance from
the ligand sphere, where the ligand sphere's radius is from the ligand center to
the outermost ligand atom. That command will produce PLANTSactiveSite.mol2 and
PLANTSactiveSiteResidues.mol2, you can use the latter as the binding site input
for PyPLIF.


Tutorial
--------

Before you start makes sure that SPORES (http://www.tcd.uni-konstanz.de/research/spores.php)
is already installed in your computer.

This tutorial is based on the ''
research (http://pubs.acs.org/doi/abs/10.1021/jm2011589)

First go to docs/ligands folder and run smitomol2.sh:

./smitomol2.sh

That script will convert the active ligands (CNS_actives_less.smi) and decoys
(Decoys_less.smi) for HRH1 (Histamine Receptor H1) to mol2 format.

The final result will be H1_ligands_ready.mol2 which will be used for retrospective
validation of PLANTS against HRH1 receptor. Next go to docs folder and run PLANTS
by entering this command:

PLANTS --mode screen plants.conf

After the docking is finished you can run pyplif by entering this command:

pyplif

OR

pyplif -c config.txt -o ifpresult.csv

The '-c' or '--config' option is for choosing the configuration file aside from
the default one (which is config.txt). While the '-o' or '--output' is for choosing
the output file, if output file has been stated in configuration file then the '-o'
option will overwrite that option.
