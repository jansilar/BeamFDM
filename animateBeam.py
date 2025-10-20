import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time

# read data
df = pd.read_csv("./build/beam_dynamic_results.csv")

# first column ... time
time_values = df.iloc[:, 0].values
# second column ... q(t)
q_values = df.iloc[:, 1].values
# third column to end ... y values of different x positions
x_labels = df.columns[2:]
# extract x positions from column labels
x = np.array([float(lbl.split('=')[1]) for lbl in x_labels])
# extract y data
Y = df.iloc[:, 2:].values  # shape: (Nt, Nx)

Nt, Nx = Y.shape
if Nt != time_values.size or Nt != q_values.size:
    raise RuntimeError("Mismatch between number of time steps and data rows")

# prepare graphs
plt.ion()  # interactive mode on
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# left graph - beam deflection (x vs y) - axes are swaped - the beam is "vertical"
line, = ax1.plot(Y[0, :], x, lw=2)
ax1.set_xlabel("y [m]")
ax1.set_ylabel("x [m]")
ax1.set_title("Beam deflection over time")
max_abs = max(abs(Y.min()), abs(Y.max()), 1e-12)
ax1.set_xlim(-1.1 * max_abs, 1.1 * max_abs)
ax1.set_ylim(x.min(), x.max())
ax1.grid(True)

# right graph: q(t)
ax2.plot(time_values, q_values, color='lightgray', lw=1)  # full curve in background
trail_line, = ax2.plot([], [], 'b-', lw=2)  # so far displayed curve
marker, = ax2.plot([], [], 'ro')  # current point
ax2.set_xlabel("Time [s]")
ax2.set_ylabel("q [N/m]")
ax2.set_title("Applied load over time")
ax2.set_xlim(time_values.min(), time_values.max())
qmin, qmax = q_values.min(), q_values.max()
qpad = max(1e-6, 0.05 * (qmax - qmin))
ax2.set_ylim(qmin - qpad, qmax + qpad)
ax2.grid(True)

# animation loop
for n, t in enumerate(time_values):
    # update left graph
    line.set_xdata(Y[n, :])

    # update right graph
    trail_line.set_data(time_values[:n+1], q_values[:n+1])
    marker.set_data([time_values[n]], [q_values[n]])

    # update titles
    ax1.set_title(f"Beam deflection — t = {t:.4f} s")
    ax2.set_title(f"Applied load — t = {t:.4f} s, q = {q_values[n]:.2f} N/m")

    # draw
    fig.canvas.draw_idle()
    plt.pause(0.04)

plt.ioff()
plt.show()