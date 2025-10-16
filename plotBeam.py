import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

print("Adresář je : " + os.getcwd())

df = pd.read_csv("./build/beam_results.csv")

xmax = np.max(np.abs(df["y[m]"]))

plt.figure(figsize=(7,4))
plt.plot(df["y[m]"]*1000, df["x[m]"], lw=2, label="Numerical")
plt.xlim(-1.1*xmax*1000, 1.1*xmax*1000)
plt.xlabel("Deflection [mm]")
plt.ylabel("height [m]")
plt.title("Beam deflection")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()