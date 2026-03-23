from Continuum_Opacity_Reader.Continuum_Opacity_Reader_test import reader as Reader

file_path = "/home/Sanhik/tardis/Dataset/s92.201.gz"  # Replace with your actual file path
reader_instance = Reader(file_path)
data = reader_instance.read()
if data is not None:
        print(data)  # Print the first few and last few rows of the DataFrame to see if the tool is working correct.