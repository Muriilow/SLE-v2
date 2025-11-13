import csv
import io
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np

# Run the likwid choosing the version of the C program and the data we want to get
def run_likwid(arg, version, option, file):

    command = "likwid-perfctr -C 0 -g " + arg + " -f -m -O ./cgSolver" 

    # Writing the output from likwid-csv in the file
    try:
        with open(file, "w") as f:
            subprocess.run(command, 
                            shell=True,
                            input=option,
                            stdout=f,
                            text=True)
    except Exception as err:
        print(f"Erro ao executar comando: {err}")

# Get the information and return only the useful ones from a dictionary 
def extract_output(file, results):

    current_region = None 
    out = ""
    with open(file, "r") as f:
        out = csv.reader(f)
        out = list(out)

        for line in out:
            if not line:
                continue
                
            # Detect Each Region
            if line[0] == 'TABLE' and 'Region' in line[1]:
                current_region = line[1].replace('Region ', '')

            elif current_region and "AVX DP MFLOP/s" in line[0] and len(line) > 1: 
                continue
            # Check if this region contains an info we might want
            elif current_region and "DP MFLOP/s" in line[0] and len(line) > 1: 
                try:
                    if current_region not in results['flops']:
                        results['flops'][current_region] = {}

                    results['flops'][current_region]['dp'] = float(line[1])
                except ValueError:
                    print("Error trying to get the content of a line")
                    pass
    
    # print(out)
    print(results)
    

# Dictionary that holds the useful information to make graphics
results = {
    'flops': {}
}
flops_dp = []
options1 = ["COM_1.in", "COM_2.in", "COM_3.in", "COM_4.in", "COM_5.in", "COM_6.in"]
options2 = ["SEM_1.in", "SEM_2.in", "SEM_3.in", "SEM_4.in", "SEM_5.in", "SEM_6.in"]


x = [32, 64, 128, 256, 512, 1000]
for file in options1:
    with open("Dados/" + file, "r") as f: 
        option = f.read()
        run_likwid("FLOPS_DP", "1", option, "tmp.txt")
        extract_output("tmp.txt", results)
        flops_dp.append(results['flops']['EXEC_1']['dp'])

plt.plot(x, flops_dp)

plt.show()
