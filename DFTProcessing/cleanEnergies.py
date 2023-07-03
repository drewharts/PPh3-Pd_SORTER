import os

#Occasionally there will be two lines of  Sum of electronic and thermal Free Energies= found
    #In this case just take the one latest in the file and adjust the energy output accordingly

def find_energies(directory):
    num_files_processed = 0
    for file_number in range(2, 500):
        filename = f"file_{file_number}.log"
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            num_files_processed += 1
            with open(file_path, 'r') as file:
                for line in file:
                    if ' Sum of electronic and thermal Free Energies=' in line:
                        energy = line.split("=")[-1].strip()
                        print(energy)
    print(f"Number of log files processed: {num_files_processed}")

def main():
    # Example usage
    directory = "/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/DFTProcessing/removedDFTSuccess"
    find_energies(directory)

if __name__ == '__main__':
    main()
