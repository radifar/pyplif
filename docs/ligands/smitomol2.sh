#!/bin/bash

babel CNS_actives_less.smi -omol2 --gen3d -p7.4 >> H1_ligands.mol2
babel Decoys_less.smi -omol2 --gen3d -p7.4 >> H1_ligands.mol2

SPORES --mode settypes H1_ligands.mol2 H1_ligands_ready.mol2
