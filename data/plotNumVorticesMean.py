import matplotlib.pyplot as plt
import numpy as np
from plotData import read_data, FONTSIZE


plt.rc("text", usetex=True)

# Access the default color cycle
color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# Get the 5th color (index 4)
fifth_color = color_cycle[4]

problem = "pointvortices/disk"
extra = "average"
quantity = "NumVorticesMeanRadius"
test = "NoRemovalVortices"


fig, ax = plt.subplots()

data = read_data("./" + problem + "/" + test + "/" + quantity + ".txt")


data[:, 0]
# normalize the y-axis by the maximum value
ax.plot(data[:, 0], data[:, 1])

data_x = data[:, 0]
data_y = 0.5 * np.sqrt(data_x)
ax.plot(
    data_x, data_y, label="$\propto\sqrt{t}$", linestyle="dashed", color=fifth_color
)

ax.legend(fontsize=FONTSIZE)
# ax.loglog()

# decorations
ax.set_xlabel("$t$", fontsize=FONTSIZE)
ax.set_ylabel(r"$\mathcal{R}_N$", fontsize=FONTSIZE)

# set ticks 0-40 to 0-pi
# ax.set_yticks([0, 5, 10, 15, 20])
# ax.set_yticklabels(
#     [
#         r"$0$",
#         r"$\pi/8$",
#         r"$\pi/4$",
#         r"$3\pi/8$",
#         r"$\pi/2$",
#     ]
# )


plt.xticks(fontsize=FONTSIZE)
plt.yticks(fontsize=FONTSIZE)


# x1, x2, y1, y2 = 0.1, 1, 0.05, 2
# ax.set_xlim(x1 - 0.04, x2 + 0.1)

plt.savefig(
    "../images/" + problem + "/" + quantity + ".pdf",
    bbox_inches="tight",
)
# plt.show()
