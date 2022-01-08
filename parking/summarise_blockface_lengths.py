import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.style.use("transport-findings.mplstyle")

if not os.path.exists("plots/"):
    os.makedirs("plots/")

with open("json/blockface_lengths.json") as f:
    params = json.load(f)

labels = ["One end", "Either end", "Middle", "Random"]

fig, ax = plt.subplots(2, 2, figsize=[6, 4])
ax = ax.flatten()

linestyles = ["-", "--", "-.", ":"]
for i,strategy in enumerate(params['strategy']):
# for i, strategy in enumerate(["one_side", "random"]):
    for j, l in enumerate(params["L"]):
        d = []
        for t in range(int(params["ergodic"] / l)):
            this_d = np.load(
                f"output/strategy_{strategy}/L_{l}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density_{t}.npy"
            )
            d.append(this_d)
        c = cm.plasma(j / len(params["L"]))
        ax[i].plot(
            np.linspace(0, params["departures_per_metre"], len(d[0])),
            np.mean(d, axis=0),
            c=c,
            rasterized=True
        )
        ax[i].fill_between(
            np.linspace(0, params["departures_per_metre"], len(d[0])),
            # y1=np.mean(d, axis=0) - np.std(d, axis=0),
            # y2=np.mean(d, axis=0) + np.std(d, axis=0),
            y1 = np.percentile(d, 20, axis=0),
            y2 = np.percentile(d, 80, axis=0),
            label=f"L={l:0.0f} m",
            fc=c,
            alpha=0.7,
            rasterized=True
        )
    ax[i].set_title(labels[i])
    ax[i].set_xscale("log")
    ax[i].set_xticks([1e-4,1e-3, 1e-2, 1e-1, 1e0, 1e1])
    ax[i].set_xlabel("Vehicles departed per metre (veh/m)")
    ax[i].set_ylabel("Linear density (-)")
    ax[i].set_xlim(xmax=params["departures_per_metre"])
    ax[i].set_ylim([0.4, 1])
    if i == 0:
        ax[i].legend(loc="lower left", framealpha=1)

# plt.subplots_adjust(left=0.09, right=0.98, bottom=0.1, top=0.95, wspace=0.25, hspace=0.45) # for 5 inch height
plt.subplots_adjust(left=0.08, right=0.98, bottom=0.11, top=0.95, wspace=0.25, hspace=0.65) # for 4 inch height
plt.savefig("plots/blockface_lengths.pdf",dpi=300)
