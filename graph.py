import subprocess
import sys
import matplotlib.pyplot as plt
import numpy as np

# Run the likwid choosing the version of the C program and the data we want to get
def run_likwid(arg, version):

    command = "likwid-perfctr -C 0 -g " + arg + " -m ./cgSolver" 
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            return result
        else:
            print(f"Erro: {result.stderr}")
            return None

    except Exception as err:
        print(f"Erro ao executar comando: {err}")
        return None

log = run_likwid("FLOPS_DP", "1")

print(log)
