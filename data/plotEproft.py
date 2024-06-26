import matplotlib.pyplot as plt
import numpy as np
from plotData import read_data, FONTSIZE


plt.rc("text", usetex=True)

problem = "embarrassed_2DNS"
kdn = "02-kdn16"
# quantity = "Energy"
quantity = "Enstrophy"
times = ["040", "080", "130", "179"]
test = "test6"
labels = ["$t = 0.4$", "$t = 0.8$", "$t = 1.3$", "$t = 1.79$"]

fig, ax = plt.subplots()

data = []
for time in times:
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

domain = 2048


for i in range(len(data)):
    data[i][:, 0] = np.arange(0, len(data[i][:, 0])) + 1
    data[i][:, 0] *= 2 * np.pi / domain
    ax.plot(data[i][:, 0], data[i][:, 1] ** 2, label=labels[i])


data_x = (np.arange(0, len(data[0][:, 0])) + 1) * (2 * np.pi / domain)
if quantity == "Energy":
    data_y = 10 / data_x**2
else:
    data_y = 2000 / data_x**2

ax.plot(data_x, data_y, label="$\propto 1/r^2$", linestyle="dashed")

ax.legend(fontsize=FONTSIZE)

ax.loglog()

# decorations
ax.set_xlabel("$r$", fontsize=FONTSIZE)
if quantity == "Energy":
    ax.set_ylabel("$E_r$", fontsize=FONTSIZE)
else:
    ax.set_ylabel("$\Omega_r$", fontsize=FONTSIZE)

# remove previous ticks
# set tick 0-pi
ax.set_xticks([np.pi / 16, np.pi / 8, np.pi / 4, np.pi / 2, np.pi])
ax.set_xticklabels([r"$\pi/16$", r"$\pi/8$", r"$\pi/4$", r"$\pi/2$", r"$\pi$"])

plt.xticks(fontsize=FONTSIZE)
plt.yticks(fontsize=FONTSIZE)

if quantity == "Energy":
    x1, x2, y1, y2 = 0.1, np.pi, 0.02, 100
    ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
    ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)
else:
    x1, x2, y1, y2 = 0.1, np.pi, 1, 100000
    ax.set_xlim(x1 - x1 / 10, x2 + x2 / 10)
    ax.set_ylim(y1 - y1 / 10, y2 + y2 / 10)

plt.savefig(
    "../images/" + problem + "/" + quantity + "_t.kdn" + kdn[-2:] + "." + test + ".pdf",
    bbox_inches="tight",
)
