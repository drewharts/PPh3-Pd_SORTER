import os

def find_first_scf_done(file_path):
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(" SCF Done: "):
                return line.strip()
    return None


def process_log_files(directory_path, output_file_path):
    with open(output_file_path, 'w') as output_file:
        for file_name in os.listdir(directory_path):
            if file_name.endswith('.log'):
                file_path = os.path.join(directory_path, file_name)
                last_scf_done = find_first_scf_done(file_path)
                if last_scf_done:
                    output_file.write(f"{file_name}: {last_scf_done}\n")

def main():
    # Provide the directory path where the log files are located
    directory_path = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemoved'

    # Provide the path for the output file
    output_file_path = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemovedEnergies.txt'

    process_log_files(directory_path, output_file_path)

if __name__ == '__main__':
        main()
