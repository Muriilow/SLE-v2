import csv
import io
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np

# Run the likwid choosing the version of the C program and the data we want to get
def run_likwid(arg, version, option, file):

    command = "likwid-perfctr -C 0 -g " + arg + " -f -m -O ./cgSolver" + version 

    # Writing the output from likwid-csv in the file
    try:
        with open(file, "w") as f:
            subprocess.run(command, 
                            shell=True,
                            input=option,
                            stdout=f,
                            text=True)
            
            print(f"{version}-------------------")
    except Exception as err:
        print(f"Erro ao executar comando: {err}")

def check_info(operation, region, results, line):
    try:
        if region not in results[operation]:
            results[operation][region] = {}

        results[operation][region] = float(line)

    except ValueError:
        print(f"Error trying to get the {operation} of a line")
        pass
# Get the information and return only the useful ones from a dictionary 
def extract_output(file, results):

    current_region = None 
    out = ""
    with open(file, "r") as f:
        out = csv.reader(f)
        out = list(out)

        isMetric = False
        for line in out:
            if not line:
                continue
                
            # Detect Each Region
            if line[0] == 'TABLE' and 'Metric' in line[2]:
                print(line)
                current_region = line[1].replace('Region ', '')
                isMetric = True
            
            elif line[0] == 'TABLE' and 'Raw' in line[2]:
                isMetric = False
                current_region = None 

            if isMetric == True and current_region and len(line) > 1:

                if "AVX DP MFLOP/s" in line[0]: 
                    check_info("flops_avx", current_region, results, line[1])

                elif "DP MFLOP/s" in line[0]: 
                    check_info("flops_dp", current_region, results, line[1])

                elif "L3 bandwidth [MBytes/s]" in line[0]: 
                    check_info("l3", current_region, results, line[1])

                elif "L2 miss ratio" in line[0]: 
                    check_info("miss", current_region, results, line[1])
    
    # print(out)
    print(results)
    

# Dictionary that holds the useful information to make graphics
results1 = {
    'flops_dp': {},
    'flops_avx': {},
    'l3': {},
    'miss': {},
}
results2 = {
    'flops_dp': {},
    'flops_avx': {},
    'l3': {},
    'miss': {},
}

# Arrays que contem cada valor para um SL de tamanho X
# O numero 1 representa a primeira versao do cgSolver 
# O numero 2 representa a segunda versao do cgSolver
flops_dp1 = []
flops_dp2 = []
flops_avx1 = []
flops_avx2 = []
mem1 = []
mem2 = []
miss1 = []
miss2 = []

# Nomes dos arquivos que possuem os argumentos de entrada
options1 = ["COM_1.in", "COM_2.in", "COM_3.in", "COM_4.in", "COM_5.in"]
# options2 = ["SEM_1.in", "SEM_2.in", "SEM_3.in", "SEM_4.in", "SEM_5.in"]

# Tamanho do SL 
x = [32, 64, 128, 256, 512]

# Para cada arquivo do array options1, vamos executar o likwid para cada tipo
# E salvar os resultados
for file in options1:
    with open("Dados/" + file, "r") as f: 
        option = f.read()

        # Running likwid with FLOPS_DP 
        run_likwid("FLOPS_DP", "1", option, "tmp.txt")
        extract_output("tmp.txt", results1)
        flops_dp1.append(results1['flops_dp']['EXEC_1'])

        run_likwid("FLOPS_DP", "2", option, "tmp.txt")
        extract_output("tmp.txt", results2)
        flops_dp2.append(results2['flops_dp']['EXEC_1'])

        # Running likwid with FLOPS_AVX 
        run_likwid("FLOPS_AVX", "1", option, "tmp.txt")
        extract_output("tmp.txt", results1)
        flops_avx1.append(results1['flops_avx']['EXEC_1'])

        run_likwid("FLOPS_AVX", "2", option, "tmp.txt")
        extract_output("tmp.txt", results2)
        flops_avx2.append(results2['flops_avx']['EXEC_1'])

        # Running likwid with L3
        run_likwid("L3", "1", option, "tmp.txt")
        extract_output("tmp.txt", results1)
        mem1.append(results1['l3']['EXEC_1'])

        run_likwid("L3", "2", option, "tmp.txt")
        extract_output("tmp.txt", results2)
        mem2.append(results2['l3']['EXEC_1'])

        # Running likwid with L2CACHE
        run_likwid("L2CACHE", "1", option, "tmp.txt")
        extract_output("tmp.txt", results1)
        miss1.append(results1['miss']['EXEC_1'])

        run_likwid("L2CACHE", "2", option, "tmp.txt")
        extract_output("tmp.txt", results2)
        miss2.append(results2['miss']['EXEC_1'])

print(flops_dp1)
print(flops_avx1)
print(mem1)
print(miss1)

plt.figure("FLOPS_DP")
plt.plot(x, flops_dp1)
plt.plot(x, flops_dp2)
plt.xscale('log')

plt.figure("FLOPS_AVX")
plt.plot(x, flops_avx1)
plt.plot(x, flops_avx2)
plt.xscale('log')

plt.figure("L3")
plt.plot(x, mem1)
plt.plot(x, mem2)
plt.xscale('log')

plt.figure("L2CACHE")
plt.plot(x, miss1)
plt.plot(x, miss2)
plt.xscale('log')
plt.show()
