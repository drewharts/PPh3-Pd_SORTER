import os

def find_number_in_logs(directory, number):
    for filename in os.listdir(directory):
        if filename.endswith(".log"):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                for line in file:
                    if str(number) in line:
                        print(f"Found in {filename}: {line.strip()}")

def main():
    # Example usage
    directory = "/Users/drewhartsfield/PycharmProjects/PPh3-Pd_SORTER/DFTProcessing/notRemovedDFTSuccess"
    number = -2647.716893
    find_number_in_logs(directory, number)

if __name__ == '__main__':
    main()
