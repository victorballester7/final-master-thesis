import matplotlib.pyplot as plt
import numpy as np
from plotData import read_data, FONTSIZE

plt.rc("text", usetex=True)

# Access the default color cycle
color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Get the 5th color (index 4)
fifth_color = color_cycle[4]

problem = "embarrassed_2DNS"
kdns = ["01-kdn8", "02-kdn16", "03-kdn32"]
quantity = "Energy"
# quantity = "Enstrophy"

test = "test7"
labels = ["$k_r = 8$", "$k_r = 16$", "$k_r = 32$"]
time = "059"  # 0.68
# domain = [2048, 2048, 4096]
domain = [4096, 4096, 4096]

fig, ax = plt.subplots()

data = []
for kdn in kdns:
    data.append(
        read_data(
            "./"
            + problem
            + "/"
            + kdn
            + "/"
            + test
            + "/"
            + quantity
            + "."
            + time
            + ".txt"
        )
    )


for i in range(len(data)):
    data[i][:, 0] = np.arange(0, len(data[i][:, 0])) + 1
    data[i][:, 0] *= 2 * np.pi / domain[i]
    ax.plot(data[i][:, 0], data[i][:, 1] ** 2, label=labels[i])

data_x = np.arange(0.01, np.pi, 0.01)
if quantity == "Energy":
    data_y = 10 / data_x**2
else:
    data_y = 2000 / data_x**2

ax.plot(data_x, data_y, label="$\propto 1/r^2$", linestyle="dashed", color=fifth_color)

ax.legend(fontsize=FONTSIZE)

ax.loglog()

# decorations
ax.set_xlabel("$r$", fontsize=FONTSIZE)
ax.set_ylabel("$E_r$", fontsize=FONTSIZE)

# remove previous ticks
# set tick 0-pi
ax.set_xticks([np.pi / 16, np.pi / 8, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
ax.set_xticklabels(
    [r"$\pi/16$", r"$\pi/8$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]
)

if quantity == "Energy":
    x1, x2, y1, y2 = 0.1, np.pi, 0.02, 100
    ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
    ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)
else:
    x1, x2, y1, y2 = 0.05, np.pi, 1, 1000000
    ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
    ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)

plt.savefig(
    "../images/" + problem + "/" + quantity + "_kdn." + test + "." + time + ".pdf",
    bbox_inches="tight",
)
