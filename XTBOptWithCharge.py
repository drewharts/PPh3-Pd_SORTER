import subprocess
import os
from decimal import Decimal

def get_total_charge(csd_code, file_path):
    with open(file_path, 'r') as file:
        found_csd_code = False
        for line in file:
            if line.startswith('CSD_code'):
                current_csd_code = line.split('=')[1].strip()
                if current_csd_code == csd_code:
                    found_csd_code = True
                else:
                    found_csd_code = False
            elif found_csd_code and line.strip() == '':
                line = previous_line  # Go to the line above the empty line
                total_charge = line.split('=')[1].strip()
                return int(Decimal(total_charge))
            previous_line = line;

    return None

def process_files(directory_path, file2_path):
    for file_name in os.listdir(directory_path):
        file_path = os.path.join(directory_path, file_name)
        if os.path.isfile(file_path):
            with open(file_path, 'r') as file1:
                print("ATTEMPTING: " + file1.name)
                csd_code = None
                for line in file1:
                    if line.startswith('CSD_code'):
                        csd_code = line.split('=')[1].strip()
                        csd_code = csd_code[:csd_code.index('|')].strip()
                        break

            if csd_code is None:
                print('CSD_code not found in the first file.')
                continue

            total_charge = get_total_charge(csd_code, file2_path)

            if total_charge is None:
                print('Total Charge not found in the second file.')
                continue

            no_ext_file = os.path.splitext(file_name)[0]
            output_file_name = f'{no_ext_file}.opt.xyz'
            output_file_path = os.path.join(directory_path, output_file_name)

            command = f'xtb {file_path} --opt --charge {total_charge}; mv xtbopt.xyz {output_file_name}'
            subprocess.run(command, shell=True,check=True)

            print(f'Processed {file_path} and saved the output as {output_file_path}')



def main():
    file1_path = "/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/smallMedLarge_removed"
    file2_path = "/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/charges.q"

    process_files(file1_path,file2_path)

if __name__ == '__main__':
        main()