#Python script to combine multiple EOS tables into one table.
#Takes files for pressure, rest-mass density, and energy density, each in the form of 2 columns:
# [pressure/rest-mass density/energy density] [enthalpy],
#with the first row being a label (i.e. the files are the output of modules_test_EoS Elliptica run).
#Writes columns to new text file with 4 columns:
#   [pressure] [rest-mass density] [energy density] [enthalpy],
#with no extra text, ready for use in Elliptica as a tabular EOS.

def combine_columns(pressure_array, rho0_array, e_array, h_array, outFileName=None):
    if (outFileName == None):
        outFileName = input("Enter name of output file: ")
        
    with open(outFileName, "w") as outFile:
        for ctr in range(len(h_array)):
            outFile.write(str(pressure_array[ctr]) + " " + str(rho0_array[ctr]) + " " + str(e_array[ctr]) + " " + str(h_array[ctr]) + "\n")
            
def read_arrays(pFileName=None, rho0FileName=None, eFileName=None):
    if (pFileName == None):
        pFileName = input("Enter name of pressure file: ")
    if (rho0FileName == None):
        rho0FileName = input("Enter name of rest-mass density file: ")
    if (eFileName == None):
        eFileName = input("Enter name of energy-density file: ")
        
    p_array = []
    rho0_array = []
    e_array = []
    h_array = []
    
    with open(pFileName, "r") as pFile:
        data = pFile.readlines()
        for point in data[1:]:
            p_array.append(point.split()[1])
            h_array.append(point.split()[0])
        
    with open(rho0FileName, "r") as rho0File:
        data = rho0File.readlines()
        for point in data[1:]:
            dataPt = point.split()
            rho0_array.append(dataPt[1])
        
    with open(eFileName, "r") as eFile:
        data = eFile.readlines()
        for point in data[1:]:
            dataPt = point.split()
            e_array.append(dataPt[1])
        
    return (p_array, rho0_array, e_array, h_array)
    
if __name__ == '__main__':
    data = read_arrays()
    combine_columns(data[0], data[1], data[2], data[3])
