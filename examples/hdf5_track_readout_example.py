import os
import glob
import matplotlib.pyplot as plt

from GAMPixTools import electron_track_tools

def main(args):
    # edep-sim generation volume is set so that z < 0, so z == 0 defines the anode
    depth = 0

    # Find files - use input args from CLI.
    track = electron_track_tools.Track(args.input_edepsim_file,
                                       args.event_index,
                                       input_format = 'dumpTree')

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
