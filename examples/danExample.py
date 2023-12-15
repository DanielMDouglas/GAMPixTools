import os
import glob
import matplotlib.pyplot as plt

from GAMPixTools import electron_track_tools

#   Some settings
# energy = 1000000
energy = 1500
# edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
depth = 0

# Find files - get list of files in folder with single energy tracks.
input_edepsim_file = '/home/dan/studies/GAMPix/muon_sample/dumpTree_single_particle_1c2374f4-aaeb-4004-983f-213dcbc63871.h5'
track = electron_track_tools.Track(input_edepsim_file,
                                   5,
                                   input_format = 'dumpTree')

# print(f'{energy/1000:3.02f} keV track {file_num:1.0f}'
#       + f', with {track.truth["num_electrons"]:4.0f} e-'
#       + f', at {depth:2.1f} m depth')

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
# print (min(track.pixel_samples['r_triggered']),
#        max(track.pixel_samples['r_triggered']))

track.display()
track.display(pixels=False)
track.display(raw_track=False)

plt.tight_layout()
plt.show()
