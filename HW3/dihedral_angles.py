import numpy as np
import mdtraj as md
#import mdshare

mol = md.load('alanine-dipeptide-nowater.pdb')
topology = mol.topology
def angle_from_positions(coordinates, indices):

    r1, r2, r3, r4 = coordinates[indices]

    u1 = r2 - r1
    u2 = r3 - r2
    u3 = r4 - r3

    angle = np.arctan2(np.linalg.norm(u2) * np.dot(u1, np.cross(u2, u3)),
                       np.dot(np.cross(u1, u2), np.cross(u2, u3)))

    return angle

def compute_dihedral_angles(molecule, i):

    C = []
    N = []
    CA = []

    for atom in molecule.topology.atoms:
        if atom.name == 'C':
            C.append(atom.index)
        elif atom.name == 'N':
            N.append(atom.index)
        if atom.name == 'CA':
            CA.append(atom.index)

    indices_phi = [C[0], N[0], CA[0], C[1]]
    indices_psi = [N[0], CA[0], C[1], N[1]]

    phi = angle_from_positions(molecule.xyz[i], indices_phi)
    psi = angle_from_positions(molecule.xyz[i], indices_psi)

    return phi, psi


phi, psi = compute_dihedral_angles(mol, 0)
print('\N{greek small letter phi}: {:.7}\n\N{greek small letter psi}: {:.7} '.format(
    phi, psi))

# Computing the dihedral angles with the mdtraj functions:

phi_mdtraj = md.compute_phi(mol)
psi_mdtraj = md.compute_psi(mol)

print('\N{greek small letter phi}: {:.7}, indices: {}\n\N{greek small letter psi}: {:.7}, indices: {}'.format(
    phi_mdtraj[1][0][0], phi_mdtraj[0][0], psi_mdtraj[1][0][0], psi_mdtraj[0][0]))

# Compute the distances of heavy atoms that forms the bonds:

distances = []
for bond in topology.bonds:
    if bond[0].element.name != 'H' and bond[1].element.name != 'H':
        d = np.linalg.norm(
            mol.xyz[0][bond[0].index] - mol.xyz[0][bond[1].index])
        distances.append(d)
        print('{}: d = {:.1} nm'.format(bond, d))
