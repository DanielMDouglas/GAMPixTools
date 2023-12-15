import os
import h5py
import numpy as np
from tqdm import tqdm 

from GAMPixTools import edepsim_tools

def main(args):
    if not os.path.isdir(args.outputPrefix):
        os.mkdir(args.outputPrefix)

    f = h5py.File(args.input)

    print ("processing events to " + args.outputPrefix)
    eventIDs = np.unique(f['segments']['eventID'])
    for thisEVID in tqdm(eventIDs):
        segmentMask = f['segments']['eventID'] == thisEVID
        evSegments = f['segments'][segmentMask]

        trajMask = f['trajectories']['eventID'] == thisEVID
        evTraj = f['trajectories'][trajMask]

        outfileName = os.path.join(args.outputPrefix, "event_"+str(thisEVID))

        edepsim_tools.quench(evSegments)

        edepsim_tools.h5_convert(evSegments, evTraj, outfileName)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input', type = str,
                        help = "input edep-sim dumped h5 to convert")
    parser.add_argument('-o', '--outputPrefix', type = str,
                        default = "tracks",
                        help = "output destination directory")
    
    args = parser.parse_args()

    main(args)
