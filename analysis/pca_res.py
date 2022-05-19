import numpy as np
from sklearn.decomposition import PCA 
from subprocess import Popen, PIPE
import subprocess
from features import invert_dist
from time import time

def get_res_coord(infilename, only_backbone=False, frame_seperator='ENDMDL'):
    #find out number of residules and number of frames
    n_residules = 0
    n_frames = 0
    atoms_per_residules = []
    res = ''
    natoms = 0
    line_cntr = 1
    with open(infilename) as infile:
        while 1:
            line = next(infile)
            w = line.split()
            if w[0] == 'ATOM':
                if line.split()[-1] == 'H' and only_backbone:
                    pass
                elif res != line[20:26]:
                    if n_residules > 0 :
                        atoms_per_residules.append(natoms)
#                        print(line_cntr, res,atoms_per_residules[-1], n_residules, len(atoms_per_residules))
                    natoms = 1
                    n_residules += 1
                    res = line[20:26]
                else:
                    natoms += 1
            #elif line.startswith('ENDMDL'):
            elif line.startswith(frame_seperator):
                atoms_per_residules.append(natoms)
#                print(n_residules,  len(atoms_per_residules))
#                print(atoms_per_residules[:10])
#                print(atoms_per_residules[-10:])
                break
            line_cntr += 1
    infile.close()
    n_frames = int(subprocess.check_output('grep '+frame_seperator+' '+infilename+' |wc -l', shell=True)[:-1])
    print(n_residules, n_frames, len(atoms_per_residules),atoms_per_residules[-5:])

    res = ''
    coords = np.zeros([n_frames, n_residules*3])
    _coords = []
    nframe = 0
    n_res = 0
    with open(infilename) as infile:
        while 1:
            try:
                line = next(infile)
                w = line[27:].split()
                if line.startswith('ATOM'):
                    if w[-1] == 'H' and only_backbone:
                        pass
                    elif res != line[20:26]:
                        if not n_res == 0 :
                            #print(line)
                            _coords = np.array(_coords)
#                            print(_coords, _coords.shape)
                            coords[nframe, (n_res-1)*3:n_res*3] = _coords.sum(axis=0)/len(_coords)
                            #print(n_res,(n_res-1)*3,n_res*3, coords[(n_res-1)*3:n_res*3,nframe], atoms_per_residules[n_res-1],_coords.shape)
    		            _coords = [[float(w[0]),float(w[1]), float(w[2])]]
                        else:
                            _coords.append([float(w[0]),float(w[1]), float(w[2])])
#                            print(n_res,n_res*3,n_res*3+2, _coords[n_res*3:n_res*3+3], atoms_per_residules[n_res])
                        n_res += 1
                	res = line[20:26]
                    else:
                        _coords.append([float(w[0]),float(w[1]), float(w[2])])
    #                    print('b',_coords_res)
    #                    _coords[natom * 3: natom * 3 + 3] = 
#                elif line.startswith('ENDMDL'):
                elif line.startswith(frame_seperator):
                    _coords = np.array(_coords)
                    coords[nframe, (n_res-1)*3:n_res*3] = _coords.sum(axis=0)/len(_coords)
                    _coords = []
                    #print(n_res,(n_res-1)*3,n_res*3, _coords[(n_res-1)*3:n_res*3], atoms_per_residules[n_res-1])
                    #print(_coords)
                    #print(len(_coords))
                    #@print(coords[:,nframe][:9],coords[:,nframe][-15:])
                    res = ''
                    if not nframe % 5: print(nframe)
                    nframe += 1
		    n_res = 0
            except StopIteration:
                print('Done getting per residule coords')
                return(coords)
    

infilename = "traj-ap.pdb"
nametoken = 'inv'
#frame_seperator = 'ENDMDL'
frame_seperator = 'END'
feature_type = 'invert_dist' # cartesian, invert_dist, compact_invert_dist
only_backbone = True
print('starting')

a = time()
print('getting residue coordinates')
coords = get_res_coord(infilename, only_backbone=only_backbone, frame_seperator=frame_seperator)
#coords = coords[:,:-9]
b = time()
print('time used for getting residual coordinates: ', b-a)
c = time()
if feature_type == 'invert_dist':
    features = invert_dist(coords)
elif feature_type == 'cartesian':
    features = coords
d = time()
print('time used for getting features: ', d-c)
#coords = coords[:-9,:]
#for _coords in coords:
#    for c in _coords:
#        print(c)
outfile1 = open('eigen_value_'+nametoken+'.log','w')
outfile2 = open('eigen_vector_top10_'+nametoken+'.log','w')
outfile3 = open('eigen_vector_top80%_'+nametoken+'.log','w')
#============================================
#outfile4 = open('mean_displacement.log','w')
#coords1 = coords.T
#coords1 = coords1.reshape([coords1.shape[0],-1,3])
#print(coords1[0,0,:])
#coords_avg = coords1.sum(axis=0) / coords1.shape[0]
#print(coords_avg.shape, coords_avg[0], coords_avg[-1])
#coords1 -= coords_avg[None,:,:]
#print(coords1[0,0,:])
#displacement = np.linalg.norm(coords1,axis=2).sum(axis=0)
#mean_displacement = displacement / coords1.shape[0]
#mean_displacement = mean_displacement / np.linalg.norm(mean_displacement)
#for v in mean_displacement:
#    outfile4.write('%10.5f\n' % (v))
#============================================
n_components = min(features.shape[0], features.shape[1])
e = time()
pca = PCA(n_components=n_components)
pca.fit(features)
f = time()
print('time used for doing PCA: ', f-e)
for ratio, value in zip(pca.explained_variance_ratio_, pca.explained_variance_):
    outfile1.write('%10.5f %10.5f \n' %(ratio, value))
acc_ratio = 0
ind = 0
for ratio in pca.explained_variance_ratio_:
    if acc_ratio < 0.8:
        acc_ratio += ratio
        ind += 1
    else:
        break


ev_sum = pca.explained_variance_ratio_.sum()
_ev_sum = 0
for i,e in enumerate(pca.explained_variance_ratio_):
    _ev_sum += e
    if _ev_sum > ev_sum * .8:
        ind_80 = i+1
print(pca.components_.shape)
if feature_type == 'invert_dist':
    pcas = pca.components_
    pca_top10 = pcas[:10]
    top10 = np.zeros([10, coords.shape[1]/3])
    top80per = np.zeros([ind_80, coords.shape[1]/3])
    for i in range(10):
    #for i in range(2):
        mat_tmp = np.zeros([coords.shape[1]/3,coords.shape[1]/3])
        ind = np.triu_indices(coords.shape[1]/3,1)
        #print(pca_top10.shape, top10.shape,mat_tmp.shape,ind[0].shape)
        mat_tmp[ind] = abs(pca_top10[i])
        mat_tmp += mat_tmp.T
        top10[i] = mat_tmp.sum(axis=0)
    for i in range(ind_80):
        mat_tmp = np.zeros([coords.shape[1]/3,coords.shape[1]/3])
        ind = np.triu_indices(coords.shape[1]/3,1)
        mat_tmp[ind] = abs(pcas[:ind_80][i])
        mat_tmp += mat_tmp.T
        top80per[i] = mat_tmp.sum(axis=0)
elif feature_type == 'cartesian':
    pcas = pca.components_.reshape([n_components,-1,3])
    top10 = np.linalg.norm(pcas[:10],axis=2)
print(pcas.shape,top10.shape)
top10_accumulated = np.dot(pca.explained_variance_ratio_[:10], top10)
top10_accumulated /= np.linalg.norm(top10_accumulated)
print(top10_accumulated.shape)
for values, value_acc in zip(top10.T, top10_accumulated):
    outfile2.write('%10.5f '*10 % (values[0],values[1],values[2],values[3],values[4],values[5],values[6],values[7],values[8],values[9]))
#    outfile2.write('%10.5f '*2 % (values[0],values[1]))
    outfile2.write('%10.5f\n' % (value_acc))

top80per_accumulated = np.dot(pca.explained_variance_ratio_[:ind_80], top80per)
top80per_accumulated /= np.linalg.norm(top80per_accumulated)
for value in top80per_accumulated:
    outfile3.write('%10.5f\n' % (value))

