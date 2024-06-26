import os
import glob
import matplotlib.pyplot as plt

from GAMPixTools import electron_track_tools

#   Some settings
# edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
depth = 0

# Find files - get list of files in folder with single energy tracks.
input_edepsim_file = '/sdf/group/neutrino/gampix/FD_size_samples/muon_1-4GeV/dumpTree/dumpTree_single_particle_5b581549-0553-4106-9266-19a6d2e49f62.h5'
track = electron_track_tools.Track(input_edepsim_file,
                                   5,
                                   input_format = 'dumpTree')

track.reset_params(charge_readout='GAMPixD')
track.readout_charge(depth)

print (track)
print (track.params.charge_readout)
print ( '  pixels charge: '
        + f'{track.pixel_samples["samples_triggered"].sum():4.0f} e-')
print ( '  coarse tiles charge: '
        + f'{track.coarse_tiles_samples["samples_triggered"].sum():4.0f} e-')
print ( '  deposited charge: '
        + f'{sum(track.raw_track["num_e"])} e-')
print (track.raw_track["num_e"].shape)
print (track.pixel_samples.keys())
print (track.pixel_samples['samples_triggered'].shape)
print (track.pixel_samples['r_triggered'].shape)
print (track.pixel_samples['samples_raw'].shape)
print (track.pixel_samples['r_raw'].shape)
print (sum(track.pixel_samples['samples_raw']))

track.display()
track.display(pixels=False)
track.display(raw_track=False)

plt.tight_layout()
plt.show()
