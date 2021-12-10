import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.style.use('transport-findings.mplstyle')

with open('json/blockface_lengths.json') as f:
    params = json.load(f)

labels = ['One end', 'Either end', 'Middle', 'Random']

fig, ax = plt.subplots(1,2,figsize=[6,2.5])
ax = ax.flatten()

linestyles = ['-','--','-.',':']
# for i,parking_mode in enumerate(params['parking_mode']):
for i,parking_mode in enumerate(['random']):
    for j,l in enumerate(params['L']):
        c = cm.plasma(j/len(params['L']))
        d = np.load(f"output/parking_mode_{parking_mode}/L_{l}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density.npy")
        ax[0].semilogx(np.linspace(0,params['departures_per_metre'],len(d)),
                    d,
                    label=f'L={l}',
                    c=c,
                    # ls=linestyles[i],
                    )

ax[0].set_xlabel('Vehicles departed per metre (veh/m)')
ax[0].set_ylabel('Linear density (-)')
# ax[0].legend(loc=0)

for i,parking_mode in enumerate(['one_side']):
    for j,l in enumerate(params['L']):
        c = cm.plasma(j/len(params['L']))
        d = np.load(f"output/parking_mode_{parking_mode}/L_{l}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/linear_density.npy")
        ax[1].semilogx(np.linspace(0,params['departures_per_metre'],len(d)),
                    d,
                    label=f'L={l:0.0f} m',
                    c=c,
                    # ls=linestyles[i],
                    )

ax[1].set_xlabel('Vehicles departed per metre (veh/m)')
ax[1].set_ylabel('Linear density (-)')
ax[1].legend(loc='center left')
# plt.savefig('plots/summary_density.pdf')

# # plt.clf()
# bins = np.linspace(3,7,20)
# bin_centres = (bins[1:] + bins[:-1])/2.
# for i,parking_mode in enumerate(parking_modes):
#     d = 2*np.load(f'output/{parking_mode}/cars.npy')
#     hist, bin_edges = np.histogram(d,
#                                    bins=bins,
#                                    # density=True,
#                                    )
#     ax[1].plot(bin_centres,
#              hist,
#              ls='-',
#              label=parking_mode.replace('_',' '),
#              )
#
# ax[1].set_xlabel('Car size (m)')
# ax[1].set_ylabel('Count (-)')
# ax[1].legend(loc=0)
# plt.savefig('plots/summary_cars.pdf')

# plt.clf()
# bins = np.logspace(-2,1,20)
# bin_centres = np.sqrt(bins[1:]*bins[:-1])
# for i,parking_mode in enumerate(params['parking_mode']):
#     d = np.load(f"output/parking_mode_{parking_mode}/L_{l}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{params['sigma_car_length'][0]}/gaps.npy")
#     hist, bin_edges = np.histogram(d,
#                                    bins=bins,
#                                    # density=True
#                                    )
#     ax[1].plot(bin_centres,
#              hist,
#              ls='-',
#              label=labels[i],
#              )
# ax[1].set_xscale('log')
# ax[1].set_xlabel('Gap size (m)')
# ax[1].set_ylabel('Count (-)')
# ax[1].legend(loc=0)

plt.subplots_adjust(left=0.09,right=0.98,bottom=0.18,top=0.99,wspace=0.25,hspace=0.3)
plt.savefig('plots/blockface_lengths.pdf')
