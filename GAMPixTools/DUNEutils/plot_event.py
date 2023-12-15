import numpy as np
import h5py
import matplotlib.pyplot as plt
import particle

mmtom = 1.e-3

def main(args):
    f = h5py.File(args.input)

    segmentMask = f['segments']['eventID'] == args.eventID
    evSegments = f['segments'][segmentMask]

    pdgs = np.unique(evSegments['pdgId'])
    colorMap = {pdgCode: plt.color_sequences['tab10'][i]
                for i, pdgCode in enumerate(pdgs)}
    labelledPDG = []
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.view_init(elev=args.elevation, azim=args.azimuth)

    for thisSegment in evSegments:
        thisPDG = thisSegment['pdgId']
        color = colorMap[thisPDG]

        if not thisPDG in labelledPDG:
            label = '$'+particle.Particle.from_pdgid(thisPDG).latex_name+'$'
            labelledPDG.append(thisPDG)
        else:
            label = None
            
        ax.plot([thisSegment['x_start']*mmtom,
                 thisSegment['x_end']*mmtom],
                [thisSegment['y_start']*mmtom,
                 thisSegment['y_end']*mmtom],
                [thisSegment['z_start']*mmtom,
                 thisSegment['z_end']*mmtom],
                color = color,
                label = label,
                )

    plt.legend(loc = 'upper right')
    ax.set_xlabel(r'x [m]')
    ax.set_ylabel(r'y [m]')
    ax.set_zlabel(r'z [m]')

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('input', type = str,
                        help = "input edep-sim dumped h5 to convert")
    parser.add_argument('-e', '--eventID', type = int,
                        default = 0,
                        help = "output destination for saved image")
    parser.add_argument('--elevation', type = float,
                        default = 30,
                        help = "view elevation angle")
    parser.add_argument('--azimuth', type = float,
                        default = 45,
                        help = "view azimuth angle")
    parser.add_argument('-o', '--output', type = str,
                        help = "output destination for saved image")
    
    args = parser.parse_args()

    main(args)
