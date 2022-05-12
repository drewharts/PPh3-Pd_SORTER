#GOALS:
    #FIGURE OUT WHAT WE WANT TO DO IN TERMS OF GOING FORWARD

    #OPTIONS:
    #1: KEEP DEVELOPING THIS TO CONSISTENTLY AND ACCURATELY REPLACE PHENYL WITH METHYL
    #2: NOW THAT I'VE LEARNED OPENBABEL BETTER, MOVE ON TO FINGERPRINTS (SHUSEN)

    #Questions to ask?
        # DO I NEED A CLASS SO THAT EACH MOLECULE IS STORED DIFFERENTLY?
        # CAN I DELETE A LIGAND THE SAME WAY I GO THROUGH A LIGAND
        # FIGURE OUT HOW TO MAKE A COPY OF THE LIGAND

from openbabel import openbabel
import time
from openbabel import pybel
import glob
import io
from PIL import ImageTk, Image
pds = 0
#TRANSITION WELL BETWEEN atomic numbers and atom names for readability
t_metals_translation = {22: "Titanium", 23: "Vanadium",24: "Chromium",25: "Manganese",26: "Iron",27: "Cobalt",28: "Nickel",29: "Copper",30: "Zinc",
                        40: "Zirconium", 41: "Niobium",42: "Molybdenum",43: "Technetium",44: "Ruthenium",45: "Rhodium",46: "Palladium",47: "Silver",48: "Cadmium",
                        72: "Hafnium", 73: "Tantalum",74: "Tungsten",75: "Rhenium",76: "Osmium",77: "Iridium",78: "Platinum",79: "Gold",80: "Mercury",
                        104: "Rutherfordium", 105: "Dubnium",106: "Seaborgium",107: "Bohrium",108: "Hassium"}
t_metals = {22: 0, 23: 0,24: 0,25: 0,26: 0,27: 0,28: 0,29: 0,30: 0,
            40: 0, 41: 0,42: 0,43: 0,44: 0,45: 0,46: 0,47: 0,48: 0,
            72: 0, 73: 0,74: 0,75: 0,76: 0,77: 0,78: 0,79: 0,80: 0,
            104: 0, 105: 0,106: 0,107: 0,108: 0}
seen_atoms = []
atoms_in_ligand = {1: 0, 6: 0, 15: 0, 0: 0}
PPh3_dict = {1: 15, 6: 18, 15: 1, 0: 0}
phenyl_dict = {1: 5, 6: 6, 15: 1, 0: 0}
atomstodelete = []
def triphenylphosphine_ligan_checker(string):
    #Counting the number of phosphoruses and phenyls
    num_p1 = string.count("c1ccccc1")
    num_p2 = string.count("c2ccccc2")
    num_p3 = string.count("c3ccccc3")
    num_P = string.count("P")

    #This is to make sure an O doesn't appear directly before or after the phosphorus
    position_p = string.find("P")
    position_o = string.find("O")

    #the number of phenyls mod 3 needs to be zero so that we know it is a triphenylphosphine
    num_p_all = num_p1 + num_p2 + num_p3
    if num_P > 0:
        if (abs(position_p - position_o) != 1):
            if (num_p_all / num_P == 3 and num_p_all % num_P == 0):
                return True

def triph_checker(string):
    length_of_string = len(string)
    num_phosphines = 0
    try:
        string_index = string.index("[P](c1ccccc1)(c1ccccc1)c1ccccc1")
        num_phosphines = string.count("[P](c1ccccc1)(c1ccccc1)(c1ccccc1)")
        num_phosphines += string.count("([P](c1ccccc1)(c1ccccc1)c1ccccc1)")
    except ValueError:
        num_phosphines = string.count( "[P](c1ccccc1)(c1ccccc1)(c1ccccc1)")
        num_phosphines += string.count("([P](c1ccccc1)(c1ccccc1)c1ccccc1)")
    else:
        string_length = len("[P](c1ccccc1)(c1ccccc1)c1ccccc1")
        if (string_index + string_length == length_of_string):
            num_phosphines +=1

    if num_phosphines > 0:
        return True

def metal_counter(molecule):
    number_of_atoms = 0
    for atom in molecule:
        number_of_atoms += 1
        if atom.OBAtom.IsMetal():
            #46 is palladium
            if atom.atomicnum == 46:
                return True
def transition_metal_counter(molecule):
    for atom in molecule:
        if atom.OBAtom.IsMetal():
            for key in t_metals:
                if key == atom.atomicnum:
                    t_metals[key] += 1
    return t_metals

def find_ligand_v2(molecule):
    ph_neighbors = 0
    b = False
    flag = False
    mol = molecule.OBMol
    print("Molecules # of atoms: ", mol.NumAtoms())
    for atom in openbabel.OBMolAtomIter(mol):
        if (atom.GetAtomicNum() == 15):
            flag = False
            seen_atoms.append(atom.GetIdx())
            atomstodelete.append(atom.GetIdx())
            atoms_in_ligand[15] += 1
            # print("Atom index: ", a.GetIdx())
            for atomChecker in openbabel.OBAtomAtomIter(atom):
                if atomChecker.GetAtomicNum() != 6 and atomChecker.GetAtomicNum() != 46:
                    flag = True
            for atomOB in openbabel.OBAtomAtomIter(atom):
                if flag:
                    atomstodelete.clear()
                    break
                # print("Bond Index of atom: ", atomOB.GetIdx())
                # if molecule.atoms[atomOB.GetIdx()].atomicnum == 6:
                if mol.GetAtom(atomOB.GetIdx()).GetAtomicNum() == 6:
                    seen_atoms.append(atomOB.GetIdx())
                    atomstodelete.append(atomOB.GetIdx())
                    atoms_in_ligand[6] += 1
                    find_ligand_v2_dfs(atomOB)
                if len(atomstodelete) == 34:
                    break

    atomstodelete.sort(reverse=True)
    atoms_to_delete = []
    for index in atomstodelete:
        mol.DeleteAtom(mol.GetAtom(index))
    print("Idxs that will be deleted",atomstodelete)
    print("Number of Idxs to be deleted: ", len(atomstodelete))
    print("Order of seen Idx:", seen_atoms)
    return mol

#and (atomOB.GetAtomicNum() == 15 or atomOB.GetAtomicNum() == 6 or atomOB.GetAtomicNum() == 1)
def find_ligand_v2_dfs(atom):
    for atomOB in openbabel.OBAtomAtomIter(atom):
        if atomOB.GetAtomicNum() == 46:
            seen_atoms.append(atomOB.GetIdx())
            atomstodelete.clear()
            break
        if atomOB.GetIdx() not in seen_atoms:
            seen_atoms.append(atomOB.GetIdx())
            atomstodelete.append(atomOB.GetIdx())
            if atomOB.GetAtomicNum() != 1 or 6 or 15:
                atoms_in_ligand[0] += 1
            else:
                atoms_in_ligand[atomOB.GetAtomicNum()] += 1
            find_ligand_v2_dfs(atomOB)
        if atomOB.GetIdx in seen_atoms:
            find_ligand_v2_dfs(atomOB)


    # print(seen_atoms)
mol_num = 0
tri_num = 0
print("Welcome")
print("What do you want to do: \n1) Sort triphenylphosphines AND output them? \n"
      "2) Sort Palladiums, get sizes of molecules, and output sorted Pds \n"
      "3) Get Number of transition metals from file\n"
      "4) Change out phenyl for methyl")
userInput = input()
if (userInput == "1"):
    start_time = time.time()
    for file in glob.glob("tmQM_X2.xyz"):
        for molecule in pybel.readfile("xyz", file):
            mol_num += 1
            smi = molecule.write(format = "smi")
            if triph_checker(smi.split()[0].strip()) == True:
                with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                    output = molecule.write("xyz")
                    f.write(output)
                # pybel.Outputfile("xyz","output.xyz",opt=None,overwrite=True).write(molecule)
                tri_num += 1
                print(smi.split()[0].strip())
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = elapsed_time // 60
    seconds = elapsed_time % 60
    print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")
    print(f"Number of Triphenylphosphines:  {tri_num}")

if (userInput == "2"):
    num_palladium_found = 0
    num_small = 0
    num_medium = 0
    num_large = 0
    num_other = 0
    start_time = time.time()
    for file in glob.glob("PPh3_2.xyz"):
        for molecule in pybel.readfile("xyz", file):
            mol_num += 1
            smi = molecule.write(format="smi")
            if (metal_counter(molecule)):

                #outputing sorted fires (if they contain Palladium from metal_counter function
                # with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                #     output = molecule.write("xyz")
                #     f.write(output)
                # pybel.Outputfile("xyz","output.xyz",opt=None,overwrite=True).write(molecule)

                #sort Palladium depending on size of molecule (small, medium, and large)
                num_palladium_found += 1
                if (molecule.OBMol.NumAtoms() > 35 and molecule.OBMol.NumAtoms() < 75):
                    num_small += 1

                elif (molecule.OBMol.NumAtoms() >= 75 and molecule.OBMol.NumAtoms() < 125):
                    num_medium += 1
                    with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                        output = molecule.write("xyz")
                        f.write(output)
                elif (molecule.OBMol.NumAtoms() >= 125):
                    num_large += 1
                else:
                    #if for some reason previous ifs didn't catch all different sizes, print molecule and how many atoms it contains
                    print(molecule.OBMol.NumAtoms())
                    num_other += 1
                    print(smi.split()[0].strip())

    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = elapsed_time // 60
    seconds = elapsed_time % 60
    print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")
    print(f"Number of Palladium:  {num_palladium_found}")
    print(f"Number of small Palladium: {num_small}")
    print(f"Number of medium Palladium: {num_medium}")
    print(f"Number of large Palladium: {num_large}")
    print(f"Number of undefined Palladium: {num_other}")

if (userInput == "3"):
    start_time = time.time()
    for file in glob.glob("PPh3-Pd_small"):
        for molecule in pybel.readfile("xyz", file):
            mol_num += 1
            smi = molecule.write(format="smi")
            results = transition_metal_counter(molecule)
    end_time = time.time()
    elapsed_time = end_time - start_time
    minutes = elapsed_time // 60
    seconds = elapsed_time % 60
    print(f"Ran through {mol_num} molecules in {minutes:.0f} minutes and {seconds:.2f} seconds.")
    for x in t_metals_translation:
        for y in t_metals:
            if x == y:
                print (t_metals_translation[x], "has", t_metals[y], "molecules")
if (userInput == "4"):
    correct_PPh3_id = 0
    incorrect_PPh3_id = 0
    one_ligand = 0
    two_ligand = 0
    start_time = time.time()
    for file in glob.glob("PPh3-Pd_medium.xyz"):
        for molecule in pybel.readfile("xyz", file):
            mol_num += 1
            print("Molecule Number:", mol_num)
            mol = pybel.Molecule(find_ligand_v2(molecule))
            if len(atomstodelete) % 34 == 0:
                if len(atomstodelete) == 34:
                    one_ligand += 1
                    with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                        output = mol.write("xyz")
                        f.write(output)
                if len(atomstodelete) == 68:
                    two_ligand += 1
                correct_PPh3_id +=1
            else:
                incorrect_PPh3_id +=1

            print('\n')

            #Atom tracker reset
            seen_atoms = []
            atomstodelete = []
            atoms_in_ligand = {1: 0, 6: 0, 15: 0, 0: 0}
    print("Ran through:", mol_num, "molecules")
    print("Correct PPh3 identification:", correct_PPh3_id)
    print("Incorrect PPh3 identification:", incorrect_PPh3_id)
    print("Correct Percentage:", format(abs(((incorrect_PPh3_id - correct_PPh3_id) / correct_PPh3_id) * 100),".2f"),"%")
    print("Number of mols with 1 PPh3:", one_ligand)
    print("Number of mosl with 2 PPh3s:", two_ligand)





