import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import hsv_to_rgb
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

plt.style.use("transport-findings.mplstyle")

# fig, ax = plt.subplots(2,2)
fig = plt.figure(figsize=[10, 6])
ax = [plt.subplot(211), plt.subplot(234), plt.subplot(235), plt.subplot(236)]

spt, spt_ax = plt.subplots(1, 1)


def make_graph(t, nt, parked_car_locs, parked_car_lengths, spot_lengths, linear_densities, L):
    ax[0].clear()
    for i in range(len(parked_car_locs)):
        ax[0].plot(
            [
                parked_car_locs[i] - parked_car_lengths[i],
                parked_car_locs[i] + parked_car_lengths[i],
            ],
            [1, 1],
            # '.-'
            "-",
        )
    ax[0].plot([0, 0], [0.95, 1.05], "k--")
    ax[0].plot([L, L], [0.95, 1.05], "k--")
    # ax[0].set_xlim([-0.1*L,1.1*L])
    ax[0].set_xlabel("Distance (m)")
    ax[0].set_ylabel("Occupancy")

    ax[1].clear()
    ax[1].semilogx(linear_densities, "k.")
    ax[1].set_xlim(xmax=nt)
    ax[1].set_xlabel("Vehicles departed")
    ax[1].set_ylabel("Linear density (-)")

    ax[2].clear()
    ax[2].hist(
        2 * np.array(parked_car_lengths),
        # bins=np.linspace(2*min_car_length,2*max_car_length,20)
        bins=20,
    )
    ax[2].set_xlabel("Car size (m)")
    ax[2].set_ylabel("Count (-)")

    ax[3].clear()
    ax[3].hist(spot_lengths, bins=np.logspace(-2, 1, 20))
    ax[3].set_xscale("log")
    ax[3].set_xlabel("Spot size (m)")
    ax[3].set_ylabel("Count (-)")

    plt.subplots_adjust(wspace=0.3, hspace=0.3, left=0.07, bottom=0.09, top=0.98, right=0.98)

    return fig


def get_random_color(id):
    rng = np.random.default_rng(seed=id)
    h = rng.random()
    s = 0.8
    v = 1.0
    return hsv_to_rgb([h, s, v])


def spatiotemporal(nt, car_ids, L):
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
    spt_ax.add_collection(collection)

    # spt_ax.set_facecolor('black') # make the axis background black

    # spt_ax.plot([parked_car_locs[i]-parked_car_lengths[i],
    #              parked_car_locs[i]+parked_car_lengths[i]],
    #              [t,t],
    #              # c=cm.plasma((car_ids[i]%10)/10),
    #              c = get_random_color(car_ids[i]["id"]),
    #              lw=0.4,
    #              # '.-'
    #              # '-'
    #              rasterized=True,
    #              )

    spt_ax.set_xlim([0, L])
    spt_ax.set_ylim([0, nt])
    spt_ax.set_xlabel("Distance (m)")
    spt_ax.set_ylabel("Vehicles departed (veh)")
    plt.subplots_adjust(bottom=0.1, left=0.1, right=0.98, top=0.98)
