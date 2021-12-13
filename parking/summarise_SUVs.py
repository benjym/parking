import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.style.use("transport-findings.mplstyle")

if not os.path.exists("plots/"):
    os.makedirs("plots/")

with open("json/SUVs.json") as f:
    params = json.load(f)

labels = ["One end", "Either end", "Middle", "Random"]

fig, ax = plt.subplots(2, 2, figsize=[6, 2.5])
ax = ax.flatten()

linestyles = ["-", "--", "-.", ":"]
# for i,strategy in enumerate(params['strategy']):
for i, strategy in enumerate(["one_side","middle","random"]):
    for j, s in enumerate(params["sigma_car_length"]):
        c = cm.plasma(j / len(params["sigma_car_length"]))
        d = np.load(
            f"output/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{s}/linear_density_0.npy"
        )
        ax[i].semilogx(
            np.linspace(0, params["departures_per_metre"], len(d)),
            d,
            label=f"$\sigma={s}$",
            c=c,
            # ls=linestyles[i],
        )

    ax[i].set_xlabel("Vehicles departed per metre (veh/m)")
    ax[i].set_ylabel("Linear density (-)")
    # ax[0].legend(loc=0)
    if i == 0: ax[i].legend(loc="center left")

bins = np.linspace(0,10,21)
bin_centres = (bins[1:] + bins[:-1])/2.
for i, strategy in enumerate(["one_side"]):
    for j, s in enumerate(params["sigma_car_length"]):
        c = cm.plasma(j / len(params["sigma_car_length"]))
        d = 2*np.load(
            f"output/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{s}/cars_0.npy"
        )
        hist, bin_edges = np.histogram(
            d,
            bins=bins,
            # density=True
        )
        ax[3].plot(
            bin_centres, hist, ls="-", label=f"$\sigma={s}$ m", c=cm.plasma(j / len(params["sigma_car_length"]))
        )

ax[3].set_xlabel("Car size (m)")
ax[3].set_ylabel("Number")
# ax[1].legend(loc='center left')

plt.subplots_adjust(left=0.09, right=0.98, bottom=0.18, top=0.99, wspace=0.25, hspace=0.3)
plt.savefig("plots/SUVs.pdf")
