import sys
import numpy as np
import matplotlib.pyplot as plt
from plotting import *

def park_the_car(parked_car_locs,test_loc,index,min_loc,max_loc):
    if parking_mode == 'test_loc':
        parked_car_locs.insert(index, test_loc)
    elif parking_mode == 'random':
        parked_car_locs.insert(index, np.random.rand()*(max_loc - min_loc) + min_loc)
    elif parking_mode == 'middle':
        parked_car_locs.insert(index, (max_loc + min_loc)/2.)
    elif parking_mode == 'one_side':
        parked_car_locs.insert(index, min_loc)
    else:
        sys.exit(f'Parking mode "{parking_mode}" is undefined')
    return parked_car_locs

def generate_new_car_length():
    # return np.random.rand()*(max_car_length - min_car_length) + min_car_length # HALF CAR LENGTH
    return np.random.normal(loc=(min_car_length+max_car_length)/2., scale=(max_car_length-min_car_length)/6.)

def get_spot_lengths(parked_car_locs,parked_car_lengths):
    if len(parked_car_locs) == 0:
        spots = [L-2*min_clearance]
    elif len(parked_car_locs) == 1:
        spots = [  parked_car_locs[0]-parked_car_lengths[0]-min_clearance,
                 L-parked_car_locs[0]-parked_car_lengths[0]-min_clearance]
    else:
        x = np.array(parked_car_locs)
        l = np.array(parked_car_lengths)
        spots = np.hstack([
            x[0] - l[0] - min_clearance,
            x[1:] - x[:-1] - l[1:] - l[:-1] - 2*min_clearance,
            L - x[-1] - l[-1]  - min_clearance
        ])
    # print(spots)
    return spots

def get_spot_bounds(index,parked_car_locs,parked_car_lengths,min_clearance,new_car_length,L):
    if len(parked_car_locs) == 0: # no cars yet
        min_loc = new_car_length + min_clearance
        max_loc = L - min_clearance - new_car_length
    elif len(parked_car_locs) == 1: # just one car
        if index == 0:
            min_loc = new_car_length + min_clearance
            max_loc = parked_car_locs[0] - parked_car_lengths[0] - 2*min_clearance - new_car_length
        else:
            min_loc = parked_car_locs[0] + parked_car_lengths[0] + 2*min_clearance + new_car_length
            max_loc = L - min_clearance - new_car_length
    else:
        if index == 0:
            min_loc = new_car_length + min_clearance
            max_loc = parked_car_locs[0] - parked_car_lengths[0] - 2*min_clearance - new_car_length
        elif index == len(parked_car_locs):
            min_loc = parked_car_locs[index-1] + parked_car_lengths[index-1] + 2*min_clearance + new_car_length
            max_loc = L - min_clearance - new_car_length
        else:
            min_loc = parked_car_locs[index-1] + parked_car_lengths[index-1] + 2*min_clearance + new_car_length
            max_loc = parked_car_locs[index]   - parked_car_lengths[index]   - 2*min_clearance - new_car_length
    return min_loc, max_loc

def is_empty(parked_car_locs,parked_car_lengths,new_car_length,test_loc):
    if len(parked_car_locs) == 0: # no cars yet
        min_loc = new_car_length + min_clearance
        max_loc = L - min_clearance - new_car_length
        parked_car_locs = park_the_car(parked_car_locs,test_loc,0,min_loc,max_loc)
        return parked_car_locs,[new_car_length],True
    elif len(parked_car_locs) == 1: # just one car
        if ( test_loc - new_car_length - min_clearance > 0 ) and ( test_loc + new_car_length + min_clearance < parked_car_locs[0] - parked_car_lengths[0] ):
            parked_car_lengths.insert(0, new_car_length)
            min_loc = new_car_length + min_clearance
            max_loc = parked_car_locs[0] - parked_car_lengths[0] - min_clearance - new_car_length
            parked_car_locs = park_the_car(parked_car_locs,test_loc,0,min_loc,max_loc)
            return parked_car_locs,parked_car_lengths,True
        elif ( test_loc - new_car_length - min_clearance > parked_car_locs[0] + parked_car_lengths[0] ) and ( test_loc + new_car_length + min_clearance < L ):
            parked_car_lengths.insert(1, new_car_length)
            min_loc = parked_car_locs[0] + parked_car_lengths[0] + min_clearance - new_car_length
            max_loc = L - min_clearance - new_car_length
            parked_car_locs = park_the_car(parked_car_locs,test_loc,1,min_loc,max_loc)
            return parked_car_locs,parked_car_lengths,True
    else:
        for i in range(len(parked_car_locs)+1):
            if i == 0:
                if ( test_loc - new_car_length - min_clearance > 0 ) and ( test_loc + new_car_length + min_clearance < parked_car_locs[i] - parked_car_lengths[i] ):
                    parked_car_lengths.insert(i, new_car_length)
                    min_loc = new_car_length + min_clearance
                    max_loc = parked_car_locs[0] - parked_car_lengths[0] - min_clearance - new_car_length
                    parked_car_locs = park_the_car(parked_car_locs,test_loc,i,min_loc,max_loc)
                    return parked_car_locs,parked_car_lengths,True
            elif i == len(parked_car_locs):
                if ( test_loc - new_car_length - min_clearance > parked_car_locs[i-1] + parked_car_lengths[i-1] ) and ( test_loc + new_car_length + min_clearance < L ):
                    parked_car_lengths.insert(i, new_car_length)
                    min_loc = parked_car_locs[i-1] + parked_car_lengths[i-1] + min_clearance + new_car_length
                    max_loc = L - min_clearance - new_car_length
                    parked_car_locs = park_the_car(parked_car_locs,test_loc,i,min_loc,max_loc)
                    return parked_car_locs,parked_car_lengths,True
            else:
                if ( test_loc - new_car_length - min_clearance > parked_car_locs[i-1] + parked_car_lengths[i-1] ) and ( test_loc + new_car_length + min_clearance < parked_car_locs[i] - parked_car_lengths[i] ):
                    parked_car_lengths.insert(i, new_car_length)
                    min_loc = parked_car_locs[i-1] + parked_car_lengths[i-1] + min_clearance + new_car_length
                    max_loc = parked_car_locs[i]   - parked_car_lengths[i]   - min_clearance - new_car_length
                    parked_car_locs = park_the_car(parked_car_locs,test_loc,i,min_loc,max_loc)
                    return parked_car_locs,parked_car_lengths,True
    # if nothing else worked
    return parked_car_locs,parked_car_lengths,False

# verbose = True
verbose = False
L = 1e4 # length of road (m)
min_car_length = 1.5 # length of HALF OF smallest car (m)
max_car_length = 3.5 # lengh of HALF OF largest car (m)
min_clearance = 0. # gap between cars at EACH END (m)
# t_start_leaving = 20

# parking_mode = 'test_loc'
# parking_mode = 'random'
# parking_mode = 'middle'
# parking_mode = 'one_side'
parking_mode = sys.argv[1]

nt = 1e5
parked_car_locs    = []
parked_car_lengths = []
linear_densities   = []
new_car_length = generate_new_car_length()
found_spot = False


for t in range(int(nt)): # remove one car every timestep

    # just remove one car
    if len(parked_car_locs) > 0:
        to_remove = np.random.randint(len(parked_car_locs))
        parked_car_locs.pop(to_remove)
        parked_car_lengths.pop(to_remove)

    # Brute force filling method
    # for i in range(1000): # aggressively try to fill every spot
    #     test_loc = L*np.random.rand()
    #     parked_car_locs,parked_car_lengths,found_spot =  is_empty(parked_car_locs,
    #                                                               parked_car_lengths,
    #                                                               new_car_length,
    #                                                               test_loc)
    #     if found_spot:
    #         new_car_length = generate_new_car_length() # another car arrives looking to park
    #         found_spot = False

    # Criteria-based filling method
    fits = True
    while fits:
        spot_lengths = get_spot_lengths(parked_car_locs,parked_car_lengths)
        avail_spots = np.array(spot_lengths) > 2*new_car_length
        num_avail_spots = np.sum(avail_spots)
        if num_avail_spots == 0: fits = False
        else:
            avail_spot_indices = np.where(avail_spots)[0]
            i = np.random.choice(avail_spot_indices)
            min_loc, max_loc = get_spot_bounds(i,parked_car_locs,parked_car_lengths,min_clearance,new_car_length,L)
            # print(max_loc-min_loc,spot_lengths[i],new_car_length)
            parked_car_locs = park_the_car(parked_car_locs,0,i,min_loc,max_loc)
            parked_car_lengths.insert(i, new_car_length)
            new_car_length = generate_new_car_length() # another car arrives looking to park —-- cars never give up, if they dont make it in this round they will persist to next round
            # print(min_loc,max_loc)
            # print(avail_spot_indices,parked_car_locs,parked_car_lengths)
    linear_densities.append( 2.*np.sum(parked_car_lengths)/L )

    if verbose:
        plt.ion()
        make_graph(t,nt,parked_car_locs,parked_car_lengths,spot_lengths,min_car_length,max_car_length,linear_densities,L)
        plt.pause(1e-3)

make_graph(t,nt,parked_car_locs,parked_car_lengths,spot_lengths,min_car_length,max_car_length,linear_densities,L)
plt.savefig(f'plots/{parking_mode}.png')

np.save(f'output/{parking_mode}_linear_density.npy',linear_densities)
np.save(f'output/{parking_mode}_spots.npy',spot_lengths)
np.save(f'output/{parking_mode}_cars.npy',parked_car_lengths)
