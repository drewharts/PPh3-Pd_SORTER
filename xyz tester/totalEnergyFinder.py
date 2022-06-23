import glob
def findTotalEnergy():
    optmol_file = glob.glob("*.opt.mol")
    with open(optmol_file) as f:
        content = f.readlines()

        thirdLine = content[3]
        print(thirdLine)
