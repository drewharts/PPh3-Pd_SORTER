import os

#This python script will iterate through the .com files from both directories and find the matching structures
#It will output the matching structures file number's in the specified output file

def compare_files(directory1, directory2, output_file_path):
    files1 = os.listdir(directory1)
    files2 = os.listdir(directory2)

    with open(output_file_path, 'w') as output_file:
        for file1 in files1:
            if file1.endswith('.com'):
                file1_path = os.path.join(directory1, file1)

                for file2 in files2:
                    if file2.endswith('.com'):
                        file2_path = os.path.join(directory2, file2)

                        if compare_lines(file1_path, file2_path):
                            output_file.write(f"{file1} - {file2}\n")

def compare_lines(file1_path, file2_path):
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2:
        lines1 = file1.readlines()
        lines2 = file2.readlines()

        for line1 in lines1:
            if line1.startswith('Pd'):
                for line2 in lines2:
                    if line2 == line1:
                        return True

    return False

def main():
    # Provide the directories to compare
    directory1 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/removed'
    directory2 = '/Users/drewhartsfield/Desktop/Ess Work/DFT/notRemoved'

    # Provide the output file path
    output_file_path = '/Users/drewhartsfield/Desktop/Ess Work/DFT/matchesFiles.txt'

    compare_files(directory1, directory2, output_file_path)

if __name__ == '__main__':
        main()
