# This is step 3 in solving problem.
#Problem this program solves:
    # Now that you have a txt file containing the matching structure's file numbers from two distinct directory...
    # You want to copy and move each file to a new directory and omit the errored DFT files

#This will use the output file created in MatchDFTfiles.py and actually grab the energies
# It will separate the energies with commas for better excel processing
# You can easily find which files didn't process correctly by uncommenting the debug line

#MAIN LOGIC - copying files from each directory into two new directoring omitting "Not Found" values

import os
import shutil

def copy_files(file_path, directory1, directory2, output_directory1, output_directory2):
    with open(file_path, 'r') as file:
        count = 1
        for line in file:
            line = line.strip()
            file_names, values = line.split('|')
            file1, file2 = file_names.split(',')

            # Skip the line if "Not Found" is present
            if 'Not Found' in values:
                continue

            file1_path = find_file(directory1, file1.strip())
            file2_path = find_file(directory2, file2.strip())

            count += 1
            if file1_path:
                copy_file(file1_path, output_directory1, count)

            if file2_path:
                copy_file(file2_path, output_directory2, count)

def find_file(directory, file_name):
    for root, dirs, files in os.walk(directory):
        if file_name in files:
            return os.path.join(root, file_name)
    return None

def copy_file(source_path, destination_directory, count):
    file_extension = os.path.splitext(source_path)[1]
    new_file_name = f"file_{count}.log"
    destination_path = os.path.join(destination_directory, new_file_name)
    shutil.copy2(source_path, destination_path)

def main():
    # Provide the file path
    file_path = '/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/DFTProcessing/matchAndPrintOutput.txt'

    # Provide the directories to search for files
    directory1 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/removed'
    directory2 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemoved'

    # Provide the output directories to copy the files
    output_directory1 = '/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/DFTProcessing/removedDFTSuccess'
    output_directory2 = '/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/DFTProcessing/notRemovedDFTSuccess'

    copy_files(file_path, directory1, directory2, output_directory1, output_directory2)

if __name__ == '__main__':
    main()
