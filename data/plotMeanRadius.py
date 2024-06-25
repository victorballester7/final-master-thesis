import matplotlib.pyplot as plt
import numpy as np
from plotData import read_data, FONTSIZE

plt.rc("text", usetex=True)

problem = "embarrassed_2DNS"
quantity = "Energy"
# quantity = "Enstrophy"

quantity += "MeanRadius"

tests = ["test5", "test6", "test7"]
kdns = ["01-kdn8", "02-kdn16", "03-kdn32"]
labels = [
    "$k_r=8$, \ \,$\mathrm{Re} = 4$",
    "$k_r=16$, $\mathrm{Re} = 8$",
    "$k_r=32$, $\mathrm{Re} = 16$",
]
domain = [1024, 2048, 4096]

fig, ax = plt.subplots()

data = []
for i in range(len(tests)):
    data.append(
        read_data(
            "./" + problem + "/" + kdns[i] + "/" + tests[i] + "/" + quantity + ".txt"
        )
    )

t_max = data[2][-1, 0]

# truncate data[0] and data[1] to match data[2]
# data[0] = data[0][data[0][:, 0] <= t_max]
# data[1] = data[1][data[1][:, 0] <= t_max]


for i in range(len(data)):
    data[i][:, 1] *= 2 * np.pi / domain[i]
    ax.plot(data[i][:, 0], data[i][:, 1] ** 2, label=labels[i])


# data_x = (np.arange(0, len(data[0][:, 0])) + 1) * (2 * np.pi / domain)
# if quantity == "Energy":
#     data_y = 10 / data_x**2
# else:
#     data_y = 2000 / data_x**2

# ax.plot(data_x, data_y, label="$\propto 1/r^2$", linestyle="dashed", color=fifth_color)

ax.legend(fontsize=FONTSIZE)

# ax.loglog()

# decorations
ax.set_xlabel("$t$", fontsize=FONTSIZE)
if quantity == "EnergyMeanRadius":
    ax.set_ylabel("$\mathcal{R}_E$", fontsize=FONTSIZE)
else:
    ax.set_ylabel("$\mathcal{R}_\Omega$", fontsize=FONTSIZE)

# remove previous ticks
# set tick 0-
if quantity == "EnergyMeanRadius":
    ax.set_yticks([0, np.pi / 16, np.pi / 8, 3 * np.pi / 16, np.pi / 4, 5 * np.pi / 16])
    ax.set_yticklabels(
        [r"$0$", r"$\pi/16$", r"$\pi/8$", r"$3\pi/16$", r"$\pi/4$", r"$5\pi/16$"],
        fontsize=FONTSIZE,
    )
else:
    ax.set_yticks([0, np.pi / 32, np.pi / 16, 3 * np.pi / 32, np.pi / 8])
    ax.set_yticklabels(
        [r"$0$", r"$\pi/32$", r"$\pi/16$", r"$3\pi/32$", r"$\pi/8$"],
        fontsize=FONTSIZE,
    )

# set font size of xticks without changing the labels
for label in ax.get_xticklabels():
    label.set_fontsize(FONTSIZE)

# if quantity == "Energy":
#     x1, x2, y1, y2 = 0.1, np.pi, 0.02, 100
#     ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
#     ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)
# else:
#     x1, x2, y1, y2 = 0.1, np.pi, 1, 100000
#     ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
#     ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)

plt.savefig(
    "../images/" + problem + "/" + quantity + ".pdf",
    bbox_inches="tight",
)
