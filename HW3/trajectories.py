import numpy as np
import mdtraj
import mdshare
import matplotlib.pyplot as plt

pdb = 'alanine-dipeptide-nowater.pdb'
files = mdshare.fetch(
    'alanine-dipeptide-*-250ns-nowater.xtc', working_directory='.')

trajs = [mdtraj.load(f, top=pdb) for f in files]

phis = [mdtraj.compute_phi(t)[1] for t in trajs]
psis = [mdtraj.compute_psi(t)[1] for t in trajs]
dihedrals = [np.hstack([_phi, _psi]) for _phi, _psi in zip(phis, psis)]
plt.hist2d(*np.concatenate(dihedrals).T, bins=100, cmap='Blues',);
plt.show()
