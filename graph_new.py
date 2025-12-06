import csv
import io
import subprocess
import os
import matplotlib.pyplot as plt
import numpy as np


# This script run two versions of cgSolver and create gaphs based on the likwid output
# To run this script you have to ensure that: 
# 1. You need to have to binaries, "cgSolver1", "cgSolver2", it need to be this exactly name
# 2. You need to have a directory caleed "Dados"
# 3. You need to have all the files exactly named as the list options1
# 4. Likwid has to installed on your machine 

# I recommend you to compile the source C code of the old version and just move to this directory

# Clear dictionary
def clear_results(results):
    for key in results:
        if isinstance(results[key], dict):
            results[key].clear()

# Create a graph with the arrays "list_old" and "list_new"
def create_graph(name, x, list_old, list_new):
    plt.figure(name)
    plt.plot(x, list_old, linestyle='--')
    plt.plot(x, list_new)
    plt.xscale('log')
    filename = os.path.join("images", f"GRAPH_{name}.png")
    plt.savefig(filename, dpi=300)


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
                current_region = line[1].replace('Region ', '')
                isMetric = True
            
            elif line[0] == 'TABLE' and 'Raw' in line[2]:
                isMetric = False
                current_region = None 

            if isMetric == True and current_region and len(line) > 1:

                if "AVX DP [MFLOP/s]" in line[0]: 
                    check_info("flops_avx", current_region, results, line[1])

                elif "DP [MFLOP/s]" in line[0]: 
                    check_info("flops_dp", current_region, results, line[1])

                elif "L3 bandwidth [MBytes/s]" in line[0]: 
                    check_info("l3", current_region, results, line[1])

                elif "L2 miss ratio" in line[0]: 
                    check_info("miss", current_region, results, line[1])
    
    if check_time == 1:
        with (
            open("timeExec1.txt", "r") as fExec,
            open("timeRes1.txt", "r") as fRes,
        ):
            linesExec = [line.strip() for line in fExec.readlines()]
            check_info("time", "EXEC_1", results, linesExec[0])

            linesRes = [line.strip() for line in fRes.readlines()]
            check_info("time", "RES_1", results, linesRes[0])

    if check_time == 2:
        with (
            open("timeExec2.txt", "r") as fExec,
            open("timeRes2.txt", "r") as fRes,
        ):
            linesExec = [line.strip() for line in fExec.readlines()]
            check_info("time", "EXEC_1", results, linesExec[0])

            linesRes = [line.strip() for line in fRes.readlines()]
            check_info("time", "RES_1", results, linesRes[0])

    # print(out)
    # print(results)

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

    flops_dp1_exec = []
    flops_dp1_res = []

    flops_dp2_exec = []
    flops_dp2_res = []

    flops_avx1_exec = []
    flops_avx1_res = []

    flops_avx2_exec = []
    flops_avx2_res = []

    mem1_exec = []
    mem1_res = []

    mem2_exec = []
    mem2_res = []

    miss1_exec = []
    miss1_res = []

    miss2_exec = []
    miss2_res = []

    time1_exec = []
    time1_res = []

    time2_exec = []
    time2_res = []

    options1 = ["COM_1.in", "COM_2.in", "COM_3.in", "COM_4.in", "COM_5.in"]

    x = [32, 64, 128, 256, 512]

    for file in options1:
        with open("Dados/" + file, "r") as f: 
            option = f.read()

            # Running likwid with FLOPS_DP 
            run_likwid("FLOPS_DP", "1", option, "tmp.txt")
            extract_output("tmp.txt", results1, False)
            flops_dp1_exec.append(results1['flops_dp']['EXEC_1'])
            flops_avx1_exec.append(results1['flops_avx']['EXEC_1'])
            flops_dp1_res.append(results1['flops_dp']['RES_1'])
            flops_avx1_res.append(results1['flops_avx']['RES_1'])

            run_likwid("FLOPS_DP", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, False)
            flops_dp2_exec.append(results2['flops_dp']['EXEC_1'])
            flops_avx2_exec.append(results2['flops_avx']['EXEC_1'])
            flops_dp2_res.append(results2['flops_dp']['RES_1'])
            flops_avx2_res.append(results2['flops_avx']['RES_1'])

            # Running likwid with L3
            run_likwid("L3", "1", option, "tmp.txt")
            extract_output("tmp.txt", results1, False)
            mem1_exec.append(results1['l3']['EXEC_1'])
            mem1_res.append(results1['l3']['RES_1'])

            run_likwid("L3", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, False)
            mem2_exec.append(results2['l3']['EXEC_1'])
            mem2_res.append(results2['l3']['RES_1'])

            # Running likwid with L2CACHE
            run_likwid("L2CACHE", "1", option, "tmp.txt")
            extract_output("tmp.txt", results1, 1)
            miss1_exec.append(results1['miss']['EXEC_1'])
            miss1_res.append(results1['miss']['RES_1'])

            time1_exec.append(results1['time']['EXEC_1'])
            time1_res.append(results1['time']['RES_1'])

            run_likwid("L2CACHE", "2", option, "tmp.txt")
            extract_output("tmp.txt", results2, 2)
            miss2_exec.append(results2['miss']['EXEC_1'])
            miss2_res.append(results2['miss']['RES_1'])

            time2_exec.append(results2['time']['EXEC_1'])
            time2_res.append(results2['time']['RES_1'])

            clear_results(results1)
            clear_results(results2)

            print(time1_exec)
            print(time1_res)
            print(time2_exec)
            print(time2_res)

    create_graph("FLOPS_DP EXECUTION", x, flops_dp1_exec, flops_dp2_exec)
    create_graph("FLOPS_DP RESIDUE", x, flops_dp1_res, flops_dp2_res)

    create_graph("FLOPS_AVX EXECUTION", x, flops_avx1_exec, flops_avx2_exec)
    create_graph("FLOPS_AVX RESIDUE", x, flops_avx1_res, flops_avx2_res)

    create_graph("L3 EXECUTION", x, mem1_exec, mem2_exec)
    create_graph("L3 RESIDUE", x, mem1_res, mem2_res)

    create_graph("L2CACHE EXECUTION", x, miss1_exec, miss2_exec)
    create_graph("L2CACHE RESIDUE", x, miss1_res, miss2_res)

    create_graph("TIME EXECUTION", x, time1_exec, time2_exec)
    create_graph("TIME RESIDUE", x, time1_res, time2_res)

    plt.close()
