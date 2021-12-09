import numpy as np
import matplotlib.pyplot as plt

# fig, ax = plt.subplots(2,2)
fig = plt.figure(figsize=[10,6])
ax = [plt.subplot(211),
      plt.subplot(234),
      plt.subplot(235),
      plt.subplot(236)
      ]


def make_graph(t,nt,parked_car_locs,parked_car_lengths,spot_lengths,min_car_length,max_car_length,linear_densities,L):
    ax[0].clear()
    for i in range(len(parked_car_locs)):
        ax[0].plot([parked_car_locs[i]-parked_car_lengths[i],
                    parked_car_locs[i]+parked_car_lengths[i]],
                    [1,1],
                    # '.-'
                    '-'
                    )
    ax[0].plot([0,0],[0.95,1.05],'k--')
    ax[0].plot([L,L],[0.95,1.05],'k--')
    # ax[0].set_xlim([-0.1*L,1.1*L])
    ax[0].set_xlabel('Distance (m)')
    ax[0].set_ylabel('Occupancy')

    ax[1].clear()
    ax[1].semilogx(linear_densities,'k.')
    ax[1].set_xlim(xmax=nt)
    ax[1].set_xlabel('Vehicles departed')
    ax[1].set_ylabel('Linear density (-)')

    ax[2].clear()
    ax[2].hist(2*np.array(parked_car_lengths),
                 bins=np.linspace(2*min_car_length,2*max_car_length,20)
                 )
    ax[2].set_xlabel('Car size (m)')
    ax[2].set_ylabel('Count (-)')

    ax[3].clear()
    ax[3].hist(spot_lengths,
                 bins=np.logspace(-2,1,20)
                 )
    ax[3].set_xscale('log')
    ax[3].set_xlabel('Spot size (m)')
    ax[3].set_ylabel('Count (-)')

    plt.subplots_adjust(wspace=0.3,hspace=0.3,left=0.07,bottom=0.09,top=0.98,right=0.98)
