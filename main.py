#UP NEXT:
    #Go on case by case basis and removing some molecules where we can't remove the ligand properly
        #Another option would be trying to remove the ligand individually on some of the molecules where the difference is > 0.1
    #After this use what we learned to do medium molecules
        #Ideally if we can get a set of about 400 good molecules
import fnmatch

from openbabel import openbabel
import time
from openbabel import pybel
import glob
import io
import os
import fnmatch
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

mol_num = 0
tri_num = 0
print("Welcome")
print("What do you want to do: \n1) Sort triphenylphosphines AND output them? \n"
      "2) Sort Palladiums, get sizes of molecules, and output sorted Pds \n"
      "3) Get Number of transition metals from file\n"
      "4) Remove PPh3\n"
      "5) Get total energy number from already calculated opt.mol files\n"
      "6) Choose this for lining up mismatched stoichs...\n"
      "7) Convert xyz files to mol files\n"
      "8) Find energy for specific stoich value")
userInput = input()
if (userInput == "1"):
    print("pls list file (with extension)")
    input = input()
    start_time = time.time()
    for file in glob.glob(input):
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
    print("pls list file (with extension)")
    input = input()
    num_palladium_found = 0
    num_small = 0
    num_medium = 0
    num_large = 0
    num_other = 0
    start_time = time.time()
    for file in glob.glob(input):
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
                    with io.open("file_" + str(mol_num) + ".mol", 'w', encoding='utf-8') as f:
                        output = molecule.write("mol")
                        f.write(output)
                elif (molecule.OBMol.NumAtoms() >= 75 and molecule.OBMol.NumAtoms() < 125):
                    num_medium += 1

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
    print("pls list file (with extension)")
    input = input()
    start_time = time.time()
    for file in glob.glob(input):
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
    print("pls list file (with extension)")
    input = input()
    correct_PPh3_id = 0
    incorrect_PPh3_id = 0
    one_ligand = 0
    two_ligand = 0
    start_time = time.time()
    for file in glob.glob(input):
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
                    with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                        output = mol.write("xyz")
                        f.write(output)
                correct_PPh3_id +=1
            else:
                incorrect_PPh3_id +=1
                with io.open("file_" + str(mol_num) + ".xyz", 'w', encoding='utf-8') as f:
                    output = mol.write("xyz")
                    f.write(output)

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

if (userInput == "5"):
    input_one = input("pls list file (containing opt.mols)")
    stoich_in = input("Do you want to print stoichs(s) or energies(e)?")
    def findTotalEnergy():
        import os
        file_counter = 0

        # Define the location of the directory
        path = r"/Users/AndrewHartsfield/PycharmProjects/WillRenameLater/" + input_one

        for filenames in os.listdir(path):
            if fnmatch.fnmatch(filenames, '*.opt.mol'):
                with open(os.path.join(path, filenames)) as myfile:
                    file_counter += 1
                    file_contents = myfile.readlines()
                    start_int = file_contents[0].find("Stoichiometry")
                    start_int_energy = file_contents[2].find("energy:")
                    end_of_stoich = 0
                    for i in range(len(file_contents[0])):
                        if i > 13:
                            if file_contents[0][i] == '|':
                                end_of_stoich = i
                    if (stoich_in == 's'):
                        print(file_contents[0][start_int + 16:end_of_stoich])
                    elif(stoich_in == 'e'):
                        print(file_contents[2][start_int_energy + 8:start_int_energy + 25])
    findTotalEnergy()

if userInput == '6':
    print("pls list file (containing opt.mols)")
    input_one = input()
    print("pls give txt file contained stoichs")
    input_two = input()
    path = r"/Users/AndrewHartsfield/PycharmProjects/WillRenameLater/" + input_one
    file = open(input_two)
    # print(file.readlines()[0])
    for i in file:
        variable = i
        for filenames in os.listdir(path):
            if fnmatch.fnmatch(filenames, '*.opt.mol'):
                with open(os.path.join(path, filenames)) as myfile:
                    file_contents = myfile.readlines()
                    start_int = file_contents[0].find("Stoichiometry")
                    start_int_energy = file_contents[2].find("energy:")
                    end_of_stoich = 0
                    for i in range(len(file_contents[0])):
                        if i > 13:
                            if file_contents[0][i] == '|':
                                end_of_stoich = i-1
                    deez = file_contents[0][start_int+16:end_of_stoich] + '\n'
                    # print(variable, deez)
                    #+ '\n'
                    if variable == deez:
                        # print(file_contents[0][start_int + 16:end_of_stoich])
                        print(file_contents[2][start_int_energy + 8:start_int_energy + 25])
                        break

if userInput == '7':
    print("pls list xyz containing all xyz's in same file")
    input = input()
    for file in glob.glob(input):
        for molecule in pybel.readfile("xyz", file):
            mol_num += 1
            with io.open("file_" + str(mol_num) + ".mol", 'w', encoding='utf-8') as f:
                output = molecule.write("mol")
                f.write(output)
    print("Run through: ", mol_num)

if userInput == '8':
    print("Please list stoich you wish to find: ")
    input_two = input()
    print("pls list file (containing opt.mols)")
    input_one = input()
    path = r"/Users/AndrewHartsfield/PycharmProjects/WillRenameLater/" + input_one

    for filenames in os.listdir(path):
        if fnmatch.fnmatch(filenames, '*.opt.mol'):
            with open(os.path.join(path, filenames)) as myfile:
                file_contents = myfile.readlines()
                start_int = file_contents[0].find("Stoichiometry")
                start_int_energy = file_contents[2].find("energy:")
                end_of_stoich = 0
                for i in range(len(file_contents[0])):
                    if i > 13:
                        if file_contents[0][i] == '|':
                            end_of_stoich = i - 1
                deez = file_contents[0][start_int + 16:end_of_stoich] + '\n'
                # print(variable, deez)
                # + '\n'
                if input_two + '\n' == deez:
                    # print(file_contents[0][start_int + 16:end_of_stoich])
                    print(filenames)
                    print(input_two, " energy is ",file_contents[2][start_int_energy + 8:start_int_energy + 25])
                    break




