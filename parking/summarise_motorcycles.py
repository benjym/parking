import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

plt.style.use("transport-findings.mplstyle")

with open("json/motorcycles.json") as f:
    params = json.load(f)

if not os.path.exists('plots/'): os.makedirs('plots/')

labels = ["One end", "Either end", "Middle", "Random"]

# for i, strategy in enumerate(params["strategy"]):
#     fig, ax = plt.subplots(1, 2, figsize=[6, 2.5])
#     ax = ax.flatten()
#     for j, m in enumerate(params["motorcycle_ratio"]):
#         d = []
#         for t in range(int(params["ergodic"] / params["L"][0])):
#             this_d = np.load(
#                 f"output/motorcycles/motorcycle_ratio_{m}/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density_{t}.npy"
#             )
#             d.append(this_d)
#         c=cm.plasma(j / len(params["motorcycle_ratio"])),
#         ax[0].fill_between(
#             np.arange(len(d[0])) / params["L"][0],
#             y1=np.mean(d, axis=0) - np.std(d, axis=0),
#             y2=np.mean(d, axis=0) + np.std(d, axis=0),
#             label=f'$m={m}$',
#             fc=c,
#             alpha=0.7,
#         )
#
#     ax[0].set_xscale('log')
#     ax[0].set_xlabel("Vehicles departed per metre (veh/m)")
#     ax[0].set_ylabel("Linear density (-)")
#     # ax[0].legend(loc=0)
#     # plt.savefig('plots/summary_density.pdf')
#
#     # plt.clf()
#     bins = np.logspace(-2, 1, 20)
#     bin_centres = np.sqrt(bins[1:] * bins[:-1])
#     for j, m in enumerate(params["motorcycle_ratio"]):
#         d = []
#         for t in range(int(params["ergodic"] / params["L"][0])):
#             this_d = np.load(
#                 f"output/motorcycles/motorcycle_ratio_{m}/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/gaps_{t}.npy"
#             )
#             d.extend(this_d)
#         # d = np.array(d).flatten()
#         hist, bin_edges = np.histogram(
#             d,
#             bins=bins,
#             # density=True
#         )
#         ax[1].plot(
#             bin_centres, hist, ls="-", label=f'$r_m={m}$', c=cm.plasma(j / len(params["motorcycle_ratio"]))
#         )
#     ax[1].set_xscale("log")
#     ax[1].set_xlabel("Gap size (m)")
#     ax[1].set_ylabel("Number")
#     ax[1].legend(loc=0)
#
#     plt.subplots_adjust(left=0.09, right=0.98, bottom=0.17, top=0.99, wspace=0.25, hspace=0.3)
#     plt.savefig(f"plots/motorcycles-{strategy}.pdf")

fig, ax = plt.subplots(2, 2, figsize=[6, 4])
ax = ax.flatten()
# for i, strategy in enumerate(["middle","random"]):
for i, strategy in enumerate(params["strategy"]):
    for j, m in enumerate(params["motorcycle_ratio"]):
        d = []
        for t in range(int(params["ergodic"] / params["L"][0])):
            this_d = np.load(
                f"output/motorcycles/motorcycle_ratio_{m}/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density_{t}.npy"
            )
            d.append(this_d)
        c=cm.plasma(j / len(params["motorcycle_ratio"])),
        ax[i].fill_between(
            np.arange(len(d[0])) / params["L"][0],
            y1=np.mean(d, axis=0) - np.std(d, axis=0),
            y2=np.mean(d, axis=0) + np.std(d, axis=0),
            label=f'$m={m}$',
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
    ax[i].set_ylim(0.7,1)
    if i == 0: ax[i].legend(loc='lower left')
    # plt.subplots_adjust(left=0.1, right=0.98, bottom=0.17, top=0.97, wspace=0.3, hspace=0.3)
    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.11, top=0.945, wspace=0.25, hspace=0.75)
    plt.savefig(f"plots/motorcycles.pdf",dpi=300)
