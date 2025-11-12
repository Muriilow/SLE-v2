import matplotlib.pyplot as plt
import numpy as np

x = np.array([2023, 2024, 2025, 2026])
y = np.array([15, 25, 30, 20])

# Creating dictionary where every arg in .plot is a key 
style = dict(marker=".", markersize="15")

plt.plot(x, y, **style)

# Customizing the Graph
plt.title("Custo", fontsize=16)
plt.xlabel("Tamanho do Sistema Linear", fontsize=12)
plt.xticks(x)
plt.grid(axis="both", linestyle="dashed")

plt.show()
