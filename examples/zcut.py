import h5py
import numpy as np
import sys

input = sys.argv[1]

with h5py.File(input, 'r+') as f:
    segs = f['segments']
    print(f'# of segs in original: {len(segs)}')
    
    mask1 = segs['z_start'][:] < 0
    mask2 = segs['z_end'][:] < 0
    mask3 = segs['z'][:] < 0
    mask = mask1 & mask2 & mask3
    cut_segs = segs[mask]

    print(f'# of segs in new: {len(cut_segs)}')
    
    del f['segments']
    f.create_dataset('segments', data=cut_segs)

    f.close()
