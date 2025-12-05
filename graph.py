import csv
import io
import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np


# This script run two versions of cgSolver and create gaphs based on the likwid output
# To run this script you have to ensure that: 
# 1. You need to have to binaries, "cgSolver1", "cgSolver2", it need to be this exactly name
# 2. You need to have a directory caleed "Dados"
# 3. You need to have all the files exactly named as the list options1
# 4. Likwid has to installed on your machine 

# I recommend you to compile the source C code of the old version and just move to this directory

# Create a graph with the arrays "list_old" and "list_new"
def create_graph(name, x, list_old, list_new):
    plt.figure(name)
    plt.plot(x, list_old, linestyle='--')
    plt.plot(x, list_new)
    plt.xscale('log')

# Run the likwid choosing the version of the C program and the data we want to get
def run_likwid(arg, version, option, file):

    command = "likwid-perfctr -C 0 -g " + arg + " -f -m -O ./cgSolver" + version 

    try:
        with open(file, "w") as f:
            subprocess.run(command, 
                            shell=True,
                            input=option,
                            stdout=f,
                            text=True)
            
        
    except Exception as err:
        print(f"Erro ao executar comando: {err}")

# Add the line content into the results dictionary 
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

    options1 = ["COM_1.in", "COM_2.in", "COM_3.in", "COM_4.in", "COM_5.in"]

    x = [32, 64, 128, 256, 512]

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
    print(time1)

    create_graph("FLOPS_DP", x, flops_dp1, flops_dp2)
    create_graph("FLOPS_AVX", x, flops_avx1, flops_avx2)
    create_graph("L3", x, mem1, mem2)
    create_graph("L2CACHE", x, miss1, miss2)
    create_graph("TEMPO", x, time1, time2)

    plt.show()
