import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

# načti data
df = pd.read_csv("./build/beam_dynamic_results.csv")

# první sloupec = čas
time_values = df.iloc[:, 0].values
# zbytek sloupců = y(x)
x_labels = df.columns[1:]

x = np.array([float(lbl.split('=')[1]) for lbl in x_labels])
Y = df.iloc[:, 1:].values  # shape: (Nt, Nx)

# příprava grafu
plt.ion()  # interaktivní mód
fig, ax = plt.subplots()
line, = ax.plot(Y[0, :], x, lw=2)
ax.set_xlabel("y [m]")
ax.set_ylabel("x [m]")
ax.set_title("Beam deflection over time")
ax.set_xlim(1.1 * Y.min(), 1.1 * Y.max())
ax.set_ylim(x.min(), x.max())

# smyčka přes časové kroky
for n, t in enumerate(time_values):
    line.set_xdata(Y[n, :])
    ax.set_title(f"t = {t:.4f} s")
    plt.pause(0.02)  # pauza ~20 ms
plt.ioff()
plt.show()