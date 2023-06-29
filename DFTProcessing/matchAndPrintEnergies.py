# This is step 2 in solving problem.
#Problem this file solves:
    #Having two directories of structures
    #One directory has full structure and another has edited structure
    #You need to keep track of which structure from one directory matches with it counterpart in the other directory.

#This python script will iterate through the txt file MatchDFTFiles outputed, it will then search for each structure's energy in its log file
#It will then output a txt file with each of the matching logs file's from each directory and their energies

#MAIN LOGIC - iterating through txt file to output energy from .com file's counterpart .log file (and outputing file names again)


import os
import shutil

def process_file(file_path, directory1, directory2):
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            file_names = line.split(' - ')
            file1 = file_names[0]
            file2 = file_names[1]

            file1_log_path = os.path.join(directory1, os.path.splitext(file1)[0] + '.log')
            file2_log_path = os.path.join(directory2, os.path.splitext(file2)[0] + '.log')

            thermal_free1 = get_thermal_free_line(file1_log_path)
            thermal_free2 = get_thermal_free_line(file2_log_path)

            print(os.path.basename(file1_log_path) + "," + os.path.basename(file2_log_path) + "|" + thermal_free1 + ", " + thermal_free2)

            # Check if the energy is found in the log files
            if thermal_free1 != 'Not Found':
                copy_file(file1_log_path, '/Users/drewhartsfield/Desktop/Ess Work/DFT/new_directory1')
            if thermal_free2 != 'Not Found':
                copy_file(file2_log_path, '/Users/drewhartsfield/Desktop/Ess Work/DFT/new_directory2')

def get_thermal_free_line(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if ' Sum of electronic and thermal Free Energies=' in line:
                return line.split("=")[-1].strip()
    return 'Not Found'

def copy_file(source_file, destination_directory):
    filename = os.path.basename(source_file)
    destination_file = os.path.join(destination_directory, filename)
    shutil.copy2(source_file, destination_file)

def main():
    # Provide the file path
    file_path = '/Users/drewhartsfield/Desktop/Ess Work/DFT/matchesFiles.txt'

    # Provide the directories to search for log files
    directory1 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/removed'
    directory2 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemoved'

    process_file(file_path, directory1, directory2)

if __name__ == '__main__':
    main()
