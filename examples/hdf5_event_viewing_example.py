import h5py
import numpy as np

# input file to examine
# should be the output of the `dumpTree` script from github.com/DUNE/larnd-sim
infile = '../neutron_sample/dumpTree_single_particle_6ba9f09c-f652-4923-a79e-45d94a346663.h5'

# open the file using h5py
f = h5py.File(infile)

# first, let's just investigate as we would with any new datafile
# at the top level, the file looks very similar to a python dictionary
print ("keys in this file:", f.keys())

# we should see three entries: 'segments', 'trajectories', and 'vertices'
# let's look at the first entry:

# this is an HDF5 dataset object, and it's not very clear what's going on
print ('segments dataset:', f['segments'])

# but we can treat it just like a numpy ndarray
# each row is a segment, and its 'columns' are defined by a custom datatype:
print ('segments datatype:', f['segments'].dtype)

# so now we encounter our first actual bits of data.  Let's look at one column:
print ('event ID\'s:', f['segments']['eventID'])

# of course, this is just an index, and there are many segments per event,
# so let's see how many individual values there are
print ('event ID\'s are just an index!', np.all(np.unique(f['segments']['eventID']) == np.arange(1000)))

# now, let's use some of these indices to look at a single event
# we'll use the usual numpy slice/mask style:
my_event_index = 123

# make a boolean mask
event_segment_mask = f['segments']['eventID'] == my_event_index
# get the segments where the mask is 'True'
event_segments = f['segments'][event_segment_mask]

# now, let's check again
print ('event ID\'s in our selection:', np.unique(event_segments['eventID'])) # should just by my_event_index

# let's do this selection for the other branches in the file

# make a boolean mask
event_trajectory_mask = f['trajectories']['eventID'] == my_event_index
# get the trajectories where the mask is 'True'
event_trajectories = f['trajectories'][event_trajectory_mask]

# trajectories contain particle-level information, like position, momentum, and time of interaction,
# as well as a 'parentID' for tracing the tree of interactions from the primary (parentID == -1)
# to the terminal particles (where energy goes below the threshold for tracking
print ('trajectory datatype:', event_trajectories.dtype)

# make a boolean mask
event_vertex_mask = f['vertices']['eventID'] == my_event_index
# get the vertices where the mask is 'True'
event_vertices = f['vertices'][event_vertex_mask]

# vertices are relatively simple for this dataset,
# but could be used to store information about primary energy, interaction type, etc
print ('vertex datatype:', event_vertices.dtype)

# let's take a look at the kinds of particles that make up this event
# we're going to look at the PDG codes for each particle
# we can get this from either 'segments' or 'trajectories', depending on
# what we want to do with the particles:
print ('unique PDG codes:', np.unique(event_segments['pdgId']))

# that's not very readable.  I like to use this package to decode them
# but it doesn't always work well for nuclei/nuclear fragments

import particle
particles = [particle.Particle.from_pdgid(this_pdg)
             for this_pdg in np.unique(event_segments['pdgId'])]
print ('particles found in this event:', [this_particle.latex_name for this_particle in particles])
