import numpy as np
import matplotlib.pyplot as plt

for parking_mode in ['one_side','middle','random']:
    d = np.load(f'output/{parking_mode}_linear_density.npy')
    plt.semilogx(d,label=parking_mode.replace('_',' '))

plt.xlabel('Vehicles departed')
plt.ylabel('Linear density (-)')
plt.legend(loc=0)
plt.savefig('plots/summary_density.pdf')

plt.clf()
for parking_mode in ['one_side','middle','random']:
    d = 2*np.load(f'output/{parking_mode}_cars.npy')
    plt.hist(d,
             bins=np.linspace(3,7,20),
             label=parking_mode.replace('_',' '),
             alpha=0.5
             )

plt.xlabel('Car size')
plt.ylabel('Count (-)')
plt.legend(loc=0)
plt.savefig('plots/summary_cars.pdf')

plt.clf()
for parking_mode in ['one_side','middle','random']:
    d = np.load(f'output/{parking_mode}_spots.npy')
    plt.hist(d,
             bins=np.logspace(-2,1,20),
             label=parking_mode.replace('_',' '),
             alpha=0.5)
plt.xscale('log')
plt.xlabel('Spot size (m)')
plt.ylabel('Count (-)')
plt.legend(loc=0)
plt.savefig('plots/summary_spots.pdf')
