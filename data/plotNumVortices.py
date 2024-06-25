import matplotlib.pyplot as plt
import numpy as np
from plotData import read_data, FONTSIZE


color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Get the 5th color (index 4)
fifth_color = color_cycle[4]

plt.rc("text", usetex=True)
plt.rc("font", family="serif")

problem = "pointvortices/disk"
extra = "average"
quantity = "NumVortices"


tests = ["R2", "R4", "R8", "R16", "R32", "R64"]
labels = ["$k_r = 8$", "$k_r=16$", "$k_r=32$", "$k_r=64$", "$k_r=128$", "$k_r=256$"]
domain = [2, 4, 8, 16, 32, 64]

fig, ax = plt.subplots()

data = []
for test in tests:
    data.append(
        read_data("./" + problem + "/" + test + "/" + extra + "/" + quantity + ".txt")
    )


for i in range(len(data)):
    # data[i][:, 0] = np.arange(0, len(data[i][:, 0])) + 1
    data[i][:, 0] /= domain[i]
    # normalize the y-axis by the maximum value
    data[i][:, 1] /= np.max(data[i][:, 1])
    if i < 4:
        ax.plot(data[i][:, 0], data[i][:, 1], label=labels[i], color=color_cycle[i])
    else:
        ax.plot(
            data[i][:, 0],
            data[i][:, 1],
            label=labels[i],
            color=color_cycle[i + 1],
        )

data_x = np.arange(0.01, 1, 0.01)
data_y = 0.2 / data_x
ax.plot(data_x, data_y, label="$\propto 1/r$", linestyle="dashed", color=fifth_color)

ax.legend(fontsize=FONTSIZE)
ax.loglog()

# decorations
ax.set_xlabel("$r$", fontsize=FONTSIZE)
ax.set_ylabel(r"$\rho_N$", fontsize=FONTSIZE)

# set ticks 0-1 to 0-pi
ax.set_xticks([0.125, 0.25, 0.5, 0.75, 1])
ax.set_xticklabels([r"$\pi/8$", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"])


plt.xticks(fontsize=FONTSIZE)
plt.yticks(fontsize=FONTSIZE)


x1, x2, y1, y2 = 0.1, 1, 0.05, 2
ax.set_xlim(x1 - 0.04, x2 + 0.1)
ax.set_ylim(y1 - 0.1, y2 + 0.1)

plt.savefig(
    "../images/" + problem + "/" + quantity + ".pdf",
    bbox_inches="tight",
)
plt.show()
