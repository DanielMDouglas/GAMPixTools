#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 20:25:38 2022

This example uses two methods of events:

    events.calculate_ordering - calculates ordering using the ckd method.  

    events.calculate_stats - calculates statistics in terms of number of events
        with full energy, missing energy, etc.   
        This method is first called in a simple form and returns simple
        scalar values
        Next it is called to generate distributions of the same values, but
        now based on the value of ckd from ordering.  This is a second
        use of calculate_stats

    We also plot the statistic distributions based on ckd values

@author: tshutt
"""

from GAMPixTools import file_tools
from GAMPixTools import events_tools
from GAMPixTools import params_tools

import os
import numpy as np
import matplotlib.pyplot as plt

#%%   Paths
paths = {}
paths['root'] = '/Users/tshutt/Documents/Work/Simulations/MEGALib'
paths['data'] = os.path.join(paths['root'], 'Data')
paths['plots'] =  os.path.join(paths['root'], 'Plots')

#%%   What to do

num_hits = 5
r_limit = 1.5  #  Upper limit on radius of first hit

num_events = 1e6
gamma_energy = 1000

source = 'FarFieldPointSource'
cos_theta = '1.0'
# cos_theta = '0.6'

geo_version = '2.0'

inc = 1

file_names = file_tools.get_file_names(
    paths,
    source,
    gamma_energy,
    cos_theta,
    geo_version,
    inc
    )

#%% Get events, apply response
events = events_tools.Events(
    file_names,
    num_events,
    read_sim_file=False,
    write_events_files=True
    )

params = params_tools.Params(file_names['path_geometry'])
events.apply_detector_response(params)

#%%     Calculate simple sum "scalar" stats and blab about it
events.calculate_stats()

print('Of ' + str(events.stats['scalar_sums_sum']) + ' total events:')
for key in events.stats['scalar_sums'].keys():
    print('   ' + key + ': ' + str(events.stats['scalar_sums'][key]))

#%%     Calculate ordering with a cut on radius,
#       then use calculate_stats to calculate stats distributions depending
#       on value of ckd

#   Make cut on events within a radius
r = np.sqrt(
    events.truth_hits['r'][0, 0, :]**2
    + events.truth_hits['r'][1, 0, :]**2
    )
r_cut = r < r_limit

#   Full cut combines r cut and selecting events with num_hits
cut = r_cut & (events.truth_hits['num_hits']==num_hits)

#   This calculates ordering, here using measured hits
events.calculate_ordering(num_hits, cut=cut)

#   Here we use stats method in events.   In this case, it includes
#   calculation of distributions based on ckd values.
#   Note that we "digitize" the ckd values and
#   feed the indices to the stats calculation.

#   First digitize ckd
num_bins = 100
ckd_edges = np.logspace(-2, 1, num_bins)
ckd_centers = ckd_edges[0:-1] + np.diff(ckd_edges) / 2
ckd_indices = np.digitize(events.ordering['ckd_value'], ckd_edges)

#   Now calculate stats distributions
events.calculate_stats(
    cut=cut,
    indices=ckd_indices,
    num_bins=ckd_centers.size
    )

sum_of_sums = np.sum(events.stats['distribution_all_sums'])

#%% Now plot

#   This next step is just done to make plotting code below tidier.
#   We create list of "sums" to be plotted, along with labels.
#   Note that the missing energy distributions include the requirement of
#   "clean entrance", which means the initial scatter is in active LAr.
#   (Events without clean entrance are not included here.)
sums_list = []
sums_list.append(
    events.stats['distribution_sums']['full_energy_ordered'] \
    + events.stats['distribution_sums']['full_energy_disordered']
    )
sums_list.append(
    events.stats['distribution_sums']['full_energy_ordered']
    )
sums_list.append(
    events.stats['distribution_sums'] \
    ['clean_entrance_missing_energy_ordered']
    )
sums_list.append(
    events.stats['distribution_sums']['full_energy_disordered']
    )
sums_list.append(
    events.stats['distribution_sums'] \
    ['clean_entrance_missing_energy_disordered']
    )

labels_list = []
labels_list.append('Full energy')
labels_list.append('Ordered, full energy')
labels_list.append('Ordered, missing energy ')
labels_list.append('Disordered, full energy')
labels_list.append('Disordered, missing energy')

fig, ax = plt.subplots()

y_max = 0.75

#  Plot cumulative distributions
for n in range(len(sums_list)):
    ax.plot(
        ckd_centers,
        np.cumsum(sums_list[n]) / sum_of_sums,
        label = labels_list[n],
        linewidth='2'
        )

ax.legend(fontsize=9)
ax.grid(linestyle=':', linewidth='1')
ax.set_xscale('log')
ax.set_yscale('linear')
ax.set_xlim(ckd_edges[0], ckd_edges[-1])
ax.set_ylim(0, y_max)
ax.set_xlabel(r'ckd value (rad$^2$)')
ax.set_ylabel(str(events.ordering['num_hits_list'][0])
                        + ' Scatter Acceptance')
title_label = f'{gamma_energy/1000:5.3f} MeV, ' \
    + str(num_hits)+ ' scatters' \
    + ', Measured hits, ' \
    + f'r<{r_limit:3.1f}'
ax.set_title(title_label)


