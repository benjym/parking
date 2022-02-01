import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def marked_spaces(L,space_length,params,ergodic): # UNTESTED!!!
    n_spaces = np.floor(L/space_length)
    spots = np.zeros(n_spaces,ergodic)
    for e in ergodic:
        for i in range(n_spaces):
            valid = False
            while not valid:
                l = np.random.normal(loc=params["mean_car_length"], scale=params["sigma_car_length"])
                if l <= space_length:
                    spots[i,e] = l
                    valid = True
    return spots


plt.style.use("transport-findings.mplstyle")

with open("json/strategies.json") as f:
    params = json.load(f)

if not os.path.exists('plots/'): os.makedirs('plots/')

labels = ["One end", "Either end", "Middle", "Random"]
fig, ax = plt.subplots(1, 2, figsize=[6, 1.75])
ax = ax.flatten()

for i, strategy in enumerate(params["strategy"]):
    d = np.load(
        f"output/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density_0.npy"
    )
    ax[0].semilogx(
        np.arange(len(d)) / params["L"],
        d,
        label=labels[i],
        c=cm.plasma(i / len(params["strategy"])),
    )

ls = ['dashed','dotted']
for i,N in enumerate([2,3]):
    # analytic solution for large L
    ax[0].semilogx(
        np.arange(len(d)) / params["L"],
        np.ones(len(d))*params['mean_car_length'][0]/(params['mean_car_length'][0] + N*params['sigma_car_length'][0]),
        ls=ls[i],
        c='k',
        # label=r'$l=\mu+'+str(N)+r'\sigma$'
    )
    # brute force method
    # space_length = params['mean_car_length'][0] + N*params['sigma_car_length'][0]
    # marked_spots = marked_spaces(L,space_length,params,ergodic)
    # ax[0].semilogx(
    # )

    ax[0].text(
        3e-1,
        0.01 + params['mean_car_length'][0]/(params['mean_car_length'][0] + N*params['sigma_car_length'][0]),
        '$l_'+str(i+1)+'=\mu+' + str(N) + r'\sigma$'
    )


# ax[0].semilogx(
#     np.arange(len(d)) / params["L"],
#     np.ones(len(d))*params['mean_car_length'][0]/(params['mean_car_length'][0] + 3*params['sigma_car_length'][0]),
#     ls='dotted',
#     c='k',
#     label=r'$\mu+3\sigma$'
# )
#
# ax[0].text(
#     1e0,
#     0.01 + params['mean_car_length'][0]/(params['mean_car_length'][0] + 3*params['sigma_car_length'][0]),
#     '$\mu+3\sigma$'
# )

ax[0].set_xlabel("Vehicles departed per metre (veh/m)")
ax[0].set_ylabel("Linear density (-)")
# ax[0].legend(loc=0)
# plt.savefig('plots/summary_density.pdf')

# plt.clf()
bins = np.logspace(-2, 1, 20)
bin_centres = np.sqrt(bins[1:] * bins[:-1])
for i, strategy in enumerate(params["strategy"]):
    d = np.load(
        f"output/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/gaps_0.npy"
    )
    hist, bin_edges = np.histogram(
        d,
        bins=bins,
        # density=True
    )
    ax[1].plot(
        bin_centres, hist, ls="-", label=labels[i], c=cm.plasma(i / len(params["strategy"]))
    )
ax[1].set_xscale("log")
ax[1].set_xlabel("Gap size (m)")
ax[1].set_ylabel("Number")
ax[1].legend(loc=0)

plt.subplots_adjust(left=0.08, right=0.99, bottom=0.24, top=0.99, wspace=0.25, hspace=0.3)
plt.savefig("plots/strategies.pdf")
