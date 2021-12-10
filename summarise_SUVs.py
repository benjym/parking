import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.style.use('transport-findings.mplstyle')

with open('json/SUVs.json') as f:
    params = json.load(f)

labels = ['One end', 'Either end', 'Middle', 'Random']

fig, ax = plt.subplots(2,2,figsize=[6,2.5])
ax = ax.flatten()

linestyles = ['-','--','-.',':']
# for i,parking_mode in enumerate(params['parking_mode']):
for i,parking_mode in enumerate(['random']):
    for j,s in enumerate(params['sigma_car_length']):
        c = cm.plasma(j/len(params['sigma_car_length']))
        d = np.load(f"output/parking_mode_{parking_mode}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{s}/linear_density.npy")
        ax[0].semilogx(np.linspace(0,params['departures_per_metre'],len(d)),
                    d,
                    label=f'$\sigma={s}$',
                    c=c,
                    # ls=linestyles[i],
                    )

ax[0].set_xlabel('Vehicles departed per metre (veh/m)')
ax[0].set_ylabel('Linear density (-)')
# ax[0].legend(loc=0)

for i,parking_mode in enumerate(['one_side']):
    for j,s in enumerate(params['sigma_car_length']):
        c = cm.plasma(j/len(params['sigma_car_length']))
        d = np.load(f"output/parking_mode_{parking_mode}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{s}/linear_density.npy")
        ax[1].semilogx(np.linspace(0,params['departures_per_metre'],len(d)),
                    d,
                    label=f'$\sigma={s}$ m',
                    c=c,
                    # ls=linestyles[i],
                    )


ax[1].set_xlabel('Vehicles departed per metre (veh/m)')
ax[1].set_ylabel('Linear density (-)')
ax[1].legend(loc='center left')

for i,parking_mode in enumerate(['one_side']):
    for j,s in enumerate(params['sigma_car_length']):
        c = cm.plasma(j/len(params['sigma_car_length']))
        d = np.load(f"output/parking_mode_{parking_mode}/L_{params['L'][0]}/mean_car_length_{params['mean_car_length'][0]}/sigma_car_length_{s}/cars.npy")
        ax[3].hist(d,
                    label=f'$\sigma={s}$ m',
                    # c=c,
                    # ls=linestyles[i],
                    )


ax[3].set_xlabel('Car size (m)')
ax[3].set_ylabel('Number (-)')
# ax[1].legend(loc='center left')

plt.subplots_adjust(left=0.09,right=0.98,bottom=0.18,top=0.99,wspace=0.25,hspace=0.3)
plt.savefig('plots/SUVs.pdf')
