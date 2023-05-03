from openbabel import openbabel
import copy
import ligandClass
import moleculeClass

def renovated_pph3_search(molecule):
    global recursive_count
    recursive_count = 0
    newMol = moleculeClass.Molecule()
    newMol.ligandArray.clear()
    ligand = ligandClass.Ligand()
    numP = 0
    mol = molecule.OBMol
    print("Molecules # of atoms: ", mol.NumAtoms())
    for atom in openbabel.OBMolAtomIter(mol):
        if atom.GetAtomicNum() == 46:
            print("Index number: " + str(atom.GetIndex()))
            for neighbor in openbabel.OBAtomAtomIter(atom):
                print("Palladium neighbor: " + str(neighbor.GetAtomicNum()))
                if neighbor.GetAtomicNum() == 15:
                    print("Index number: " + str(neighbor.GetIndex()))
                    numP += 1
                    recursive_phosphine(neighbor,newMol,ligand)
    print("Number of Phospherus connected to Palladium: " + str(numP))
    if (len(newMol.ligandArray) == 2):
        newLigand = newMol.ligandArray[0]
        print(len(newLigand.atomObjects))
        for atom in newLigand.atomObjects:
            mol.DeleteAtom(atom)
        return mol,newLigand
    elif (len(newMol.ligandArray) == 1):
        newLigand = newMol.ligandArray[0]
        print(len(newLigand.atomObjects))
        for atom in newLigand.atomObjects:
            mol.DeleteAtom(atom)
        return mol
    return mol

seenPPh3AtomsObjects = []
seenPPh3AtomsNum = []
def recursive_phosphine(atom,newMol,ligand):
    global recursive_count
    if recursive_count == 0:
        seenPPh3AtomsObjects.clear()
        seenPPh3AtomsNum.clear()
    recursive_count += 1
    if atom not in seenPPh3AtomsObjects:
        for atomOB in openbabel.OBAtomAtomIter(atom):
            print("Index number: " + str(atomOB.GetIndex()))
            seenPPh3AtomsObjects.append(atomOB)
            seenPPh3AtomsNum.append(atomOB.GetAtomicNum())
            recursive_phosphine(atomOB,newMol,ligand)
    #first need to check if we've seen all the neighbors or not
    print(seenPPh3AtomsNum)
    print("num 6 atoms: " + str(seenPPh3AtomsNum.count(6)))
    print("num 46 atoms: " + str(seenPPh3AtomsNum.count(46)))
    if (seenPPh3AtomsNum.count(6) == 3 and seenPPh3AtomsNum.count(46) == 1):
        for atom in seenPPh3AtomsObjects:
            if atom.GetAtomicNum() != 46:
                recurse_p_group(atom)
        if (seenPPh3AtomsNum.count(1) == 15 and seenPPh3AtomsNum.count(6) == 18 and seenPPh3AtomsNum.count(46) == 1 and seenPPh3AtomsNum.count(15) == 1):
            #you need to remove Pd so that it doesn't get included in the ligand information
            for atom in seenPPh3AtomsObjects:
                if atom.GetAtomicNum() == 46:
                    seenPPh3AtomsObjects.remove(atom)
            ligand.atomObjects = copy.copy(seenPPh3AtomsObjects)
            newMol.ligandArray.append(ligand)
            seenPPh3AtomsObjects.clear()
            seenPPh3AtomsNum.clear()
        else:
            seenPPh3AtomsNum.clear()
            seenPPh3AtomsObjects.clear()



def recurse_p_group(atom):
    for atomOB in openbabel.OBAtomAtomIter(atom):
        if atomOB not in seenPPh3AtomsObjects:
            print("Index number: " + str(atomOB.GetIndex()))
            seenPPh3AtomsObjects.append(atomOB)
            seenPPh3AtomsNum.append(atomOB.GetAtomicNum())
            recurse_p_group(atomOB)
    print("seen atoms: " + str(seenPPh3AtomsNum))
