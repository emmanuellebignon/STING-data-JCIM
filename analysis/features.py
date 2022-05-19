import numpy as np
from numpy.linalg import norm
from time import time

def invert_dist(coords):
    #shape of the coords is nframes by natoms*3
    #shape of the features is nframes by (natoms * (natoms - 1) / 2)

    #feats = np.empty((self.nfeatures))
    #idx = 0
    #for n1, coords1 in enumerate(conf):
    #    for n2 in range(n1 + 1, self.natoms):
    #        coords2 = conf[n2]
    #        feats[idx] = 1 / np.linalg.norm(coords1 - coords2 + self._delta)
    #        idx += 1
    print('===========')
    coords = coords.reshape([coords.shape[0], coords.shape[1]/3,3])
    nframes = coords.shape[0]
    natoms = coords.shape[1]
    features = np.zeros([nframes, natoms * (natoms - 1) / 2])
    for i in range(nframes):
        ind_b = 0
        ind_e = natoms - 1
        for j in range(natoms-1):
            tmp = coords[i,j] - coords[i,j+1:natoms]
            try:
                _feats = 1 / norm(coords[i,j] - coords[i,j+1:natoms],axis=1)
            except RuntimeWarning:
                print('zero division',i,j)
                quit()

#            print(ind_b, ind_e, ind_e - ind_b, len(_feats))
#            print(features[i, ind_b:ind_e].shape, _feats.shape, features.shape)
            features[i, ind_b:ind_e] = _feats
            ind_b = ind_e
            ind_e += natoms - j - 2
        if not i % 50: print(i)

    return features




#def compact_invert_dist()
    #shape of the coords is nframes by natoms*3
    #shape of the features is nframes by (4 * (natoms - 4) + 6
