import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

# načti data
df = pd.read_csv("./build/beam_dynamic_results.csv")

# první sloupec = čas
time_values = df.iloc[:, 0].values
# druhý sloupec = q (externí síla)
q_values = df.iloc[:, 1].values
# od třetího sloupce dál = y(x) pro různé x
x_labels = df.columns[2:]
x = np.array([float(lbl.split('=')[1]) for lbl in x_labels])
Y = df.iloc[:, 2:].values  # shape: (Nt, Nx)

Nt, Nx = Y.shape
if Nt != time_values.size or Nt != q_values.size:
    raise RuntimeError("Mismatch between number of time steps and data rows")

# příprava grafu
plt.ion()  # interaktivní mód
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# levý graf: průhyb prutu (y vs x)
line, = ax1.plot(Y[0, :], x, lw=2)
ax1.set_xlabel("y [m]")
ax1.set_ylabel("x [m]")
ax1.set_title("Beam deflection over time")
max_abs = max(abs(Y.min()), abs(Y.max()), 1e-12)
ax1.set_xlim(-1.1 * max_abs, 1.1 * max_abs)
ax1.set_ylim(x.min(), x.max())
ax1.grid(True)

# pravý graf: q(t) — zobrazíme celou závislost šedě a postupně doplňujeme bod/stopu
ax2.plot(time_values, q_values, color='lightgray', lw=1)  # úplná křivka jako pozadí
trail_line, = ax2.plot([], [], 'b-', lw=2)  # dosud vykreslené body q
marker, = ax2.plot([], [], 'ro')  # aktuální bod
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("q [N/m]")
ax2.set_title("Applied load over time")
ax2.set_xlim(time_values.min(), time_values.max())
qmin, qmax = q_values.min(), q_values.max()
qpad = max(1e-6, 0.05 * (qmax - qmin))
ax2.set_ylim(qmin - qpad, qmax + qpad)
ax2.grid(True)

# animační smyčka
for n, t in enumerate(time_values):
    # update levý graf
    line.set_xdata(Y[n, :])

    # update pravý graf (trail + marker)
    trail_line.set_data(time_values[:n+1], q_values[:n+1])
    marker.set_data([time_values[n]], [q_values[n]])

    # aktualizovat titulky
    ax1.set_title(f"Beam deflection — t = {t:.4f} s")
    ax2.set_title(f"Applied load — t = {t:.4f} s, q = {q_values[n]:.2f} N/m")

    # vykreslit
    fig.canvas.draw_idle()
    plt.pause(0.04)

plt.ioff()
plt.show()