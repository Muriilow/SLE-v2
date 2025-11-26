import csv
import io
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np

def create_graph(name, x, list_old, list_new):
    plt.figure(name)
    plt.plot(x, list_old)
    plt.plot(x, list_new)
    plt.xscale('log')

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
def extract_output(file, results, check_time):
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
    
    if check_time == True:
        with open("time.txt", "r") as f:
            lines = [line.strip() for line in f.readlines()]
            results["time"] = float(lines[0])

    # print(out)
    print(results)

if __name__ == "__main__":
    # Dictionary that holds the useful information to make graphics
    results1 = {
        'flops_dp': {},
        'flops_avx': {},
        'l3': {},
        'miss': {},
        'time': {}
    }
    results2 = {
        'flops_dp': {},
        'flops_avx': {},
        'l3': {},
        'miss': {},
        'time': {}
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
    time1 = []
    time2 = []

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
            extract_output("tmp.txt", results1, False)
            flops_dp1.append(results1['flops_dp']['EXEC_1'])
            flops_avx1.append(results1['flops_avx']['EXEC_1'])

            run_likwid("FLOPS_DP", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, False)
            flops_dp2.append(results2['flops_dp']['EXEC_1'])
            flops_avx2.append(results2['flops_avx']['EXEC_1'])

            # Running likwid with L3
            run_likwid("L3", "1", option, "tmp.txt")
            extract_output("tmp.txt", results1, False)
            mem1.append(results1['l3']['EXEC_1'])

            run_likwid("L3", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, False)
            mem2.append(results2['l3']['EXEC_1'])

            # Running likwid with L2CACHE
            run_likwid("L2CACHE", "1", option, "tmp.txt")
            extract_output("tmp.txt", results1, True)
            miss1.append(results1['miss']['EXEC_1'])
            time1.append(results1['time'])

            run_likwid("L2CACHE", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, True)
            miss2.append(results2['miss']['EXEC_1'])
            time2.append(results2['time'])


    print(flops_dp1)
    print(flops_avx1)
    print(mem1)
    print(miss1)

    create_graph("FLOPS_DP", x, flops_dp1, flops_dp2)
    create_graph("FLOPS_AVX", x, flops_avx1, flops_avx2)
    create_graph("L3", x, mem1, mem2)
    create_graph("L2CACHE", x, miss1, miss2)
    create_graph("TEMPO", x, time1, time2)

    plt.show()
