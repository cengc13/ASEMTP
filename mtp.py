"""ASE calculator interface for Moment Tensor Potential."""

import os
# import time
import numpy as np
# import subprocess

from ase.calculators.calculator import Calculator, all_changes
from utils import atoms_to_cfg, read_cfg

class MTP(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    MTP_EXE = 'mlp'

    nolabel = True

    default_parameters = {'unique_elements': ['Fe', 'Co', 'Ni', 'Cr', 'Al'],
                          'tmp_folder': '/tmp/'}


    def __init__(self, mtp='pot.mtp', **kwargs):
        self.mtp = mtp
        Calculator.__init__(self, mtp, **kwargs)

    def initialize(self, atoms):
        self.numbers = atoms.get_atomic_numbers()

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if 'numbers' in system_changes:
            self.initialize(self.atoms)

        unique_elements = self.parameters.unique_elements

        file = os.path.join(self.parameters.tmp_folder, 'in.cfg')
        output = os.path.join(self.parameters.tmp_folder, 'out.cfg')
        # file = 'in.cfg' # for test
        # output = 'out.cfg' # for test
        _atoms = atoms.copy()
        atoms_to_cfg(_atoms, file, unique_elements)

        os.system(f'mlp calc-efs {self.mtp} {file} {output}')
        # p=subprocess.Popen([self.MTP_EXE, 'calc-efs', self.mtp, file, output],
        #                  stdout=subprocess.PIPE) # faster than above
        # time.sleep(0.1)
        energy, forces, stress = read_cfg(output, unique_elements)
        self.results['energy'] = energy
        self.results['forces'] = np.array(forces)
        self.results['stress'] = np.array(stress)

