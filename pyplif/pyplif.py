#!/usr/bin/env python

import getopt, sys
from openbabel import OBMol, OBConversion
from time import time
from tanimoto_coef import *
from ring import *
from interactions import *


if __name__ == "__main__":
    x = time()
    

    #  Default configuration
    config    = "config.txt"
    dl        = "cetirizine.mol2_entry_00001_conf_01.mol2"
    dp        = "cetirizine.mol2_entry_00001_conf_01_protein.mol2"
    pbf       = "protein_bindingsite_fixed.mol2"
    pr        = "H1_site.mol2"
    lr        = "Doxepin.mol2"
    multiconf = "no"


    options, args = getopt.getopt(sys.argv[1:], '', ['dl=', 'dp=', 'pbf=',
                    'config=', 'multiconf'])

    for option, value in options:
        if option == '--dl':
            dl = value
        elif option == '--dp':
            dp = value
        elif option == '--pbf':
            pbf = value
        elif option == '--config':
            config = value
        elif option == '--multiconf':
            multiconf = value

    residuechoice = ['ASP107', 'TRP158', 'PHE199', 'TRP428', 'PHE432', 'PHE435']
    

    # opening the molecule files
    conv1, conv2, conv3 = OBConversion(), OBConversion(), OBConversion()
    conv1.SetInFormat("mol2")
    conv2.SetInFormat("mol2")
    conv3.SetInFormat("mol2")

    protfix  = OBMol()
    dockprot = OBMol()
    docklig  = OBMol()
    protref  = OBMol()
    ligref   = OBMol()

    conv1.ReadFile(protfix, pbf)
    conv1.ReadFile(protref, pr)
    conv1.ReadFile(ligref, lr)
    dockprottemp = conv2.ReadFile(dockprot, dp)
    dockligtemp  = conv3.ReadFile(docklig,  dl)

    refresdict    = getresiduedict(protref)
    fixresdict    = getresiduedict(protfix)
    refringdict   = getringdict(protref)
    fixringdict   = getringdict(protfix)

    a = time()
    print "Time needed for preparation is %.3f s." % (a-x)
    ringinteraction(refresdict, refringdict, residuechoice, protref, ligref)
    ringinteraction(fixresdict, fixringdict, residuechoice, protfix, docklig)
    
    b = time()
    print "Time needed for identifying pi-interactions is %.3f s." % (b-a)
    
    otherinteractions(refresdict, residuechoice, protref, ligref)
    otherinteractions(fixresdict, residuechoice, protfix, docklig)
    hbonddockprot(fixresdict, residuechoice, dockprot, docklig)

    c = time()
    print "Time needed for identifying other interactions is %.3f s." % (c-b)


    if multiconf == "yes":
        dockprotconf = []
        dockligconf  = []
        while dockligtemp:
            dockprotconf.append(dockprot)
            dockligconf.append(docklig)
            dockprot = OBMol()
            docklig  = OBMol()
            dockprottemp = conv2.Read(dockprot)
            dockligtemp  = conv3.Read(docklig)

        confnum = 1
        for confp, confl in zip(dockprotconf, dockligconf):
            print "test"
    
    tc = gettcfromdict(refresdict, fixresdict, residuechoice)
    
    print "        %s %s" % (ligref.GetTitle()[:10], docklig.GetTitle()[:10])
    for residue in residuechoice:
        print residue, " %s  %s" % (str(refresdict[residue])[10:-2], str(fixresdict[residue])[10:-2])
    print "Tanimoto Coefficient: %.3f" % tc
                    
    y = time()
    print 'Total time taken %.3f s.' % (y-x)





