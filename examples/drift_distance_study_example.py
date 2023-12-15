import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import tqdm

from GAMPixTools import electron_track_tools

# edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
depth = 0 # include no additional shifting

# Find files - get list of files in folder with single energy tracks.
input_edepsim_file = '/home/dan/studies/GAMPix/muon_sample/dumpTree_single_particle_1c2374f4-aaeb-4004-983f-213dcbc63871.h5'

zbins = np.linspace(-10, 0, 51)
zbinCenters = 0.5*(zbins[1:] + zbins[:-1])

zMeasurements = []
chargeRatioMeasurements = []
chargeRatioMeasurements_err = []

# evID = 6
# for evID in range(2,10):
for evID in tqdm.tqdm(range(100)):
    track = electron_track_tools.Track(input_edepsim_file,
                                       evID,
                                       input_format = 'dumpTree')
    
    track.reset_params(charge_readout='GAMPixD')
    track.params.charge_drift['drift_length'] = 0.000000001
    print(track.params.charge_drift['drift_length'])
    track.readout_charge(depth)

    # print (track)
    # print (track.params.charge_readout)
    # print ( '  pixels charge: '
    #         + f'{track.pixel_samples["samples_triggered"].sum():4.0f} e-')
    # print ( '  coarse tiles charge: '
    #         + f'{track.coarse_tiles_samples["samples_triggered"].sum():4.0f} e-')
    # print ( '  deposited charge: '
    #         + f'{sum(track.raw_track["num_e"])} e-')
    # print (np.min(track.pixel_samples['r_triggered'][2,:]))
    # print (np.max(track.pixel_samples['r_triggered'][2,:]))
    # print (np.min(track.raw_track['r'][2,:]))
    # print (np.max(track.raw_track['r'][2,:]))

    raw_counts, bins = np.histogram(track.raw_track['r'][2,:],
                                    weights = track.raw_track['num_e'],
                                    bins = zbins)

    # print (hist)

    # trig_counts, bins = np.histogram(track.pixel_samples['r_triggered'][2,:],
    #                                  weights = track.pixel_samples['samples_triggered'],
    #                                  bins = zbins)
    trig_counts, bins = np.histogram(track.pixel_samples['r_raw'][2,:],
                                     weights = track.pixel_samples['samples_raw'],
                                     bins = zbins)

    # print (hist)

    mask = raw_counts != 0.
    # zMeasurements = np.concatenate((zMeasurements,
    #                                 zbinCenters[mask]))
    # chargeRatioMeasurements = np.concatenate((chargeRatioMeasurements,
    #                                           trig_counts[mask]/raw_counts[mask]))
    zMeasurements = np.concatenate((zMeasurements,
                                    np.array([np.mean(track.raw_track['r'][2,:])])))
    chargeRatioMeasurements = np.concatenate((chargeRatioMeasurements,
                                              np.array([np.sum(track.drifted_track['num_e'])/np.sum(track.raw_track['num_e'])])))
    # print (trig_counts[mask]/raw_counts[mask])
    # # print (zMeasurements)
    # if np.any(trig_counts[mask]/raw_counts[mask] > 1):
    #     print ('something weird here')

    #     track.display()
    #     track.display(pixels=False)
    #     track.display(raw_track=False)

    #     plt.show()
        
plt.scatter(zMeasurements, chargeRatioMeasurements)
plt.show()
