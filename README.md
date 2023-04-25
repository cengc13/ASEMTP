# ASEMTP

ASE calculator for [moment tensor potential](https://gitlab.com/ashapeev/mlip-2).

## Install

Need to install [ase](https://gitlab.com/ase/ase) and [mlip-2](https://gitlab.com/ashapeev/mlip-2). Current MTP version is under development at Skoltech and led by [Alex Shapeev](http://www.shapeev.com). Please contact MTP authors for the software and make citations as appropriate.

## Usage

```python
from mtp import MTP
from ase.io import read

calc = MTP(mtp='pot.mtp', unique_elements=['Fe', 'Co', 'Ni', 'Cr', 'Al'])

atoms = read('example.traj')
atoms = atoms.repeat(6)
atoms.calc = calc
# print(len(atoms))

print(atoms.get_forces())
print(atoms.get_potential_energy())
print(atoms.get_stress())
# atoms.write('sp.traj')
```