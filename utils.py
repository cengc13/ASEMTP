import re
import numpy as np
from ase.cell import Cell

def atoms_to_cfg(atoms, file, unique_elements):
    write_f, write_e = True, True
    f = open(file, 'w')
    ele_dict = {ele.capitalize(): int(i) for i, ele in
                enumerate(unique_elements)}
    try:
        e = atoms.get_potential_energy()
    except:
        write_e = False
    f.write('BEGIN_CFG\n')
    f.write(' Size\n')
    size = len(atoms)
    f.write(f'    {int(size)}\n')
    f.write(' Supercell\n')
    cell = atoms.get_cell()
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[0][0], cell[0][1], cell[0][2]))
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[1][0], cell[1][1], cell[1][2]))
    f.write("{0:<9}{1}      {2}      {3}\n".format(
        '', cell[2][0], cell[2][1], cell[2][2]))
    try:
        fs = atoms.get_forces()
    except:
        write_f = False
    if write_f:
        f.write(' AtomData:  id type       cartes_x      cartes_y'
                '      cartes_z           fx          fy          fz\n')
    else:
        f.write(' AtomData:  id type       cartes_x      cartes_y'
                '      cartes_z\n')
    pos = atoms.positions
    symbols = atoms.symbols
    for i in range(size):
        aid = int(i+1)
        atype = ele_dict[symbols[i]]
        x, y, z = pos[i]
        if write_f:
            f_x, f_y, f_z = fs[i]
            f.write('{0:>14}{1:>5}{2:>16.8f}{3:>16.8f}{4:>16.8f}{5:>12.6f}{6:>12.6f}{7:>12.6f}\n'.format(
                aid, atype, x, y, z, f_x, f_y, f_z))
        else:
            f.write('{0:>14}{1:>5}{2:>16.8f}{3:>16.8f}{4:>16.8f}\n'.format(aid, atype, x, y, z))
    if write_e:
        f.write(' Energy\n')
        f.write(f'{e:16.6f}\n')
    f.write('END_CFG\n')
    f.write('\n')
    f.close()
    return

def read_cfg(file, symbols):
    """Heavily adapted from `mlearn` package:
    https://github.com/materialsvirtuallab/mlearn"""
    with open(file, 'r') as f:
        lines = f.read()
    block_pattern = re.compile('BEGIN_CFG\n(.*?)\nEND_CFG', re.S)
    # size_pattern = re.compile('Size\n(.*?)\n SuperCell', re.S | re.I)
    lattice_pattern = re.compile('SuperCell\n(.*?)\n AtomData', re.S | re.I)
    position_pattern = re.compile('fz\n(.*?)\n Energy', re.S)
    energy_pattern = re.compile('Energy\n(.*?)\n (?=PlusStress|Stress)', re.S)
    stress_pattern = re.compile('xy\n(.*?)(?=\n|$)', re.S)
    formatify = lambda string: [float(s) for s in string.split()]
    for block in block_pattern.findall(lines):
        # size_str = size_pattern.findall(block)[0]
        # size = int(size_str.lstrip())
        lattice_str = lattice_pattern.findall(block)[0]
        lattice = np.array(list(map(formatify, lattice_str.split('\n'))))
        cell = Cell(lattice)
        volume = cell.volume
        position_str = position_pattern.findall(block)[0]
        position = np.array(list(map(formatify, position_str.split('\n'))))
        # types = position[:, 1].tolist()
        # chemical_symbols = [symbols[int(ind)] for ind in types]
        forces = position[:, 5:8].tolist()
        position = position[:, 2:5]
        energy_str = energy_pattern.findall(block)[0]
        energy = float(energy_str.lstrip())
        stress_str = stress_pattern.findall(block)[0]
        virial_stress = np.array(list(map(formatify, stress_str.split()))).reshape(6, )/volume
        virial_stress = virial_stress.tolist()
    return energy, forces, virial_stress
