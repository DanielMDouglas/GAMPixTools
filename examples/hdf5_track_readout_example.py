import os
import glob
import matplotlib.pyplot as plt

from GAMPixTools import electron_track_tools

def main(args):
    # edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
    depth = 0

    # Find files - get list of files in folder with single energy tracks.
    # input_edepsim_file = '/home/dan/studies/GAMPix/muon_sample/dumpTree_single_particle_1c2374f4-aaeb-4004-983f-213dcbc63871.h5'
    # input_edepsim_file = '/home/dan/studies/gpt_test/neutron_sample/dumpTree_single_particle_6ba9f09c-f652-4923-a79e-45d94a346663.h5'
    track = electron_track_tools.Track(args.input_edepsim_file,
                                       5,
                                       input_format = 'dumpTree')

    # print(f'{energy/1000:3.02f} keV track {file_num:1.0f}'
    #       + f', with {track.truth["num_electrons"]:4.0f} e-'
    #       + f', at {depth:2.1f} m depth')

    # set the readout to GAMPix for DUNE
    track.reset_params(charge_readout='GAMPixD')
    track.readout_charge(depth)

    print ("track object", track)
    print ("charge readout measurements", track.params.charge_readout)
    print ( '  pixels charge: '
            + f'{track.pixel_samples["samples_triggered"].sum():4.0f} e-')
    print ( '  coarse tiles charge: '
            + f'{track.coarse_tiles_samples["samples_triggered"].sum():4.0f} e-')
    print ( '  deposited charge: '
            + f'{sum(track.raw_track["num_e"])} e-')
    print ("number of charge bundles simulated:", track.raw_track["num_e"].shape)
    print ("number of electrons per bundle:", track.raw_track["num_e"])
    print ("pixel sample keys:", track.pixel_samples.keys())
    print ("number of triggered pixels:", track.pixel_samples['samples_triggered'].shape)
    print ("triggered pixel position measurements:", track.pixel_samples['r_triggered'])
    print ("charge arrived at anode:", track.pixel_samples['samples_raw'])
    print ("charge positions at anode:", track.pixel_samples['r_raw'])
    print ("total charge at anode:", sum(track.pixel_samples['samples_raw']))
    # print (min(track.pixel_samples['r_triggered']),
    #        max(track.pixel_samples['r_triggered']))

    track.display()
    track.display(pixels=False)
    track.display(raw_track=False)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input_edepsim_file',
                        type = str,
                        help = 'input file from which to read and simulate an event')
    parser.add_argument('-e', '--event_index',
                        type = int,
                        default = 5,
                        help = 'index of the event within the input file to be simulated')

    args = parser.parse_args()

    main(args)
