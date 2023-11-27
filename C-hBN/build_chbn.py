from ase.io import read, write
import numpy as np


def Set_group_from_bonds(atom):
    i = neighbor_list('i', atom, 3)
    coord = np.bincount(i)
    b2_ids = coord <= 2
    group = list(map(int, b2_ids))
    atom.set_array('group', group, int, ())
    return atom


# df: fixed length; db: broast length
# ng: number of group without db, dy
# hdir: direction of heat flux.  0:x, 1:y, 2:z
def Set_element(atom, df=10, db=10, ng=10, hdir=0, reps=[1, 1, 1]):
    cent_mass = atom.get_center_of_mass()
    symbols = atom.get_chemical_symbols()
    pos = atom.get_positions()[:, hdir]
    nat = atom.get_global_number_of_atoms()
    for n in range(nat):
        if pos[n] > cent_mass[hdir]+df*0.5+0.5:
            symbols[n] = 'C'
    atom.set_chemical_symbols(symbols)
    atom.set_pbc([1, 1, 0])

    out_atom = atom.repeat(reps)
    cells = out_atom.get_cell()[:3, :3]
    pos = out_atom.get_positions()[:, hdir]
    nat = out_atom.get_global_number_of_atoms()

    group = np.zeros((nat))
    dg = (cells[hdir][hdir]-df-db*2) / ng
    for n in range(nat):
        if pos[n] < df:
            gi = 0
        elif pos[n] < df+db:
            gi = 1
        elif pos[n] > cells[hdir][hdir]-db:
            gi = 2
        else:
            gi = (pos[n]-df-db) // dg + 3
        if df == 0: gi -= 1
        if db == 0: gi -= 2
        if df == 0 and db == 0: gi = (pos[n]-dg/2+0.4) // dg
        if gi < 0:  gi += ng
        if gi > ng: gi -= ng
        group[n] = gi
    out_atom.set_array('group', group, int, ())

    return out_atom



for repx in [2, 4]:
    rep = [46*2, 21, 1]
    ngs = int(rep[0]/23) * repx
    edg = 'zig'
    unit_cell = read(f'POSCAR.{edg}')
    cbn_rep = unit_cell.repeat(rep)
    out_mod = Set_element(cbn_rep, df=0, db=0, ng=ngs, hdir=0, reps=[repx, 1, 1])
    out_mod.write(f'ChBN-Rx{rep[0]}-Ry{rep[1]}-Ng{ngs}-{edg}-Lat{repx}.xyz', format='extxyz')

for repx in [2]:
    rep = [46*4, 21, 1]
    ngs = int(rep[0]/23) * repx
    edg = 'zig'
    unit_cell = read(f'POSCAR.{edg}')
    cbn_rep = unit_cell.repeat(rep)
    out_mod = Set_element(cbn_rep, df=0, db=0, ng=ngs, hdir=0, reps=[repx, 1, 1])
    out_mod.write(f'ChBN-Rx{rep[0]}-Ry{rep[1]}-Ng{ngs}-{edg}-Lat{repx}.xyz', format='extxyz')

exit()

for ix in [160, 320, 640]:
    rep = [ix, 12, 1]
    ngs = int(rep[0]/40)  # 80 -> 4
    edg = 'arm'
    unit_cell = read(f'POSCAR.{edg}')
    cbn_rep = unit_cell.repeat(rep)
    Set_element(cbn_rep, df=0, db=0, ng=ngs, hdir=0, reps=[repx, 1, 1])
    cbn_rep.write(f'ChBN-Rx{rep[0]}-Ry{rep[1]}-Ng{ngs}-{edg}.xyz', format='extxyz')


# Rx80-Ry12-Ng4-arm.xyz 3840
# Rx160-Ry12-Ng8-arm.xyz 7680
# Rx320-Ry12-Ng16-arm.xyz 15360
# Rx640-Ry12-Ng32-arm.xyz 30720

# Rx46-Ry21-Ng4-zig.xyz 3864
# Rx92-Ry21-Ng8-zig.xyz 7728
# Rx184-Ry21-Ng16-zig.xyz 15456
# Rx368-Ry21-Ng32-zig.xyz 30912
