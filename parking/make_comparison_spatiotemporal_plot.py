import sys, pickle, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import hsv_to_rgb
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
# from plotting import spatiotemporal

plt.style.use("transport-findings.mplstyle")

def spatiotemporal(ax, nt, car_ids, L, params):
    patches = []
    for i in range(len(car_ids)):
        # [bottom left], width, height
        if "t_f" in car_ids[i]:
            rect = mpatches.Rectangle(
                [car_ids[i]["loc"] - car_ids[i]["len"], car_ids[i]["t_i"]],
                2 * car_ids[i]["len"],
                car_ids[i]["t_f"] - car_ids[i]["t_i"],
                ec="none",
            )
        else:
            rect = mpatches.Rectangle(
                [car_ids[i]["loc"] - car_ids[i]["len"], car_ids[i]["t_i"]],
                2 * car_ids[i]["len"],
                nt - car_ids[i]["t_i"],
                ec="none",
            )
        patches.append(rect)

    colors = np.random.rand(len(patches))
    collection = PatchCollection(patches,
                                 # cmap = plt.cm.plasma,
                                 cmap = plt.cm.rainbow,
                                 # cmap = plt.cm.hsv,
                                 # alpha=0.3,
                                 )
    collection.set_array(colors)
    ax.add_collection(collection)

    # spt_ax.set_facecolor('black') # make the axis background black

    ax.set_xlim([0, L])
    ax.set_ylim([0, nt])

names = ['One end', 'Either end', 'Middle', 'Random']

spt, spt_ax = plt.subplots(1, 4, figsize = [6,3])

json_file = sys.argv[1]
with open(json_file) as f:
    params = json.load(f)

nt = int(np.ceil(params["departures_per_metre"] * params["L"][0]))
for i,strategy in enumerate(params["strategy"]):

    params["outfolder"] = f"output/strategy_{strategy}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}"

    print(params["outfolder"])
    with open(f"{params['outfolder']}/car_ids.pickle", "rb") as f:
        car_ids = pickle.load(f)

    spatiotemporal(spt_ax[i], nt, car_ids, params["L"][0], params)

    spt_ax[i].set_title(names[i])
    if i > 0:
        spt_ax[i].set_xticks([])
        spt_ax[i].set_yticks([])

spt_ax[0].set_xlabel("Distance (m)")
spt_ax[0].set_ylabel("Vehicles departed (veh)")

plt.subplots_adjust(bottom=0.14, left=0.085, right=0.995, top=0.93)
plt.savefig('plots/spatiotemporal_all.pdf')
