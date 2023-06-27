import os
#This will use the output file created in MatchDFTfiles.py and actually grab the energies
# It will separate the energies with commas for better excel processing
# You can easily find which files didn't process correctly by uncommenting the debug line

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

            print(thermal_free1 + ", " + thermal_free2)
            #for debugging purposes
            # print(file1+ ":" + thermal_free1 + " " + file2+ ":" + thermal_free2)

def get_thermal_free_line(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if ' Sum of electronic and thermal Free Energies=' in line:
                return line.split("=")[-1].strip()
    return 'Not Found'

def main():
    # Provide the file path
    file_path = '/Users/drewhartsfield/Desktop/Ess Work/DFT/matchesFiles.txt'

    # Provide the directories to search for log files
    directory1 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/removed'
    directory2 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemoved'

    process_file(file_path, directory1, directory2)

if __name__ == '__main__':
        main()
