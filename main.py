#GOALS:
    #Possibly refactor required to create class that can easily store a dictionary or some sort of data structure
    #to count the number of each transitional metal.

    #yes we are focusing on Pd but there are possibly more common ones usable out of other 2000 structures

    #Two plausible outcomes from investigating...
        # 1: my triphenylphosphine algorithm is too picky and is not pulling enough triphenylphosphines out of database
        # 2: The database just doesn't have that many triphenylsophine ligands attached to Palladium

    #Questions to ask?
        # How many triphenylphosphine ligands attached to Palladium do we expect?

from openbabel import openbabel
import time
from openbabel import pybel
import glob
import io
pds = 0
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
    num_phosphines = string.count( "[P](c1ccccc1)(c1ccccc1)c1ccccc1")

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
mol_num = 0
tri_num = 0
print("Welcome")
print("What do you want to do: \n1) Sort triphenylphosphines AND output them? \n2) Sort Palladiums, get sizes of molecules, and output sorted Pds")
userInput = input()
if (userInput == "1"):
    start_time = time.time()
    for file in glob.glob("fileName.xyz"):
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
    for file in glob.glob("fileName.xyz"):
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





