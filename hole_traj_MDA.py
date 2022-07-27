#!/usr/bin/python
import MDAnalysis.transformations.translate
import os
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import hole2
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def chain_traj_xtract(resid_range, chain, universe):
    prot = universe.select_atoms(resid_range)
    prot.write(f"./chain{chain}/{chain}coord.pdb")
    with MDAnalysis.Writer(f"./chain{chain}/{chain}traj.xtc", prot.n_atoms) as W:
        for ts in universe.trajectory:
            W.write(prot)


def chain_rms_fit(chain, ref_res_list):
    corrRefRes = [res - 1 for res in ref_res_list]
    raw_traj = f"./chain{chain}/{chain}traj.xtc"
    coord = f"./chain{chain}/{chain}coord.pdb"
    ref_frame = mda.Universe(coord, coord)
    raw_univ = mda.Universe(coord, raw_traj)
    alignment = align.AlignTraj(
        raw_univ, ref_frame, filename=f"./chain{chain}/{chain}rmsfit.dcd"
    )
    alignment.run()
    align_traj = f"./chain{chain}/{chain}rmsfit.dcd"
    align_univ = mda.Universe(coord, align_traj)
    refAtoms = align_univ.residues[corrRefRes].atoms
    center_transformation = MDAnalysis.transformations.center_in_box(
        refAtoms, center="mass", point=[0, 0, 0]
    )
    align_univ.trajectory.add_transformations(center_transformation)
    prot = align_univ.select_atoms("all")
    with MDAnalysis.Writer(f"./chain{chain}/{chain}centered.dcd", prot.n_atoms) as W:
        for ts in align_univ.trajectory:
            W.write(prot)
    return align_univ


def run_hole_traj(chain, ref_res, step=100):
    fit = chain_rms_fit(chain, ref_res)
    trj_len = len(fit.trajectory)
    print("Traj len", trj_len)
    H = hole2.HoleAnalysis(
        fit,
        cpoint=[0, 0, 0],
        cvect=[0, 0, 1],
        sample=0.125,
        executable="~/hole2/exe/hole",
        end_radius=4.0,
        prefix=f"./chain{chain}/",
    )
    H.run(step=step)
    profile = H.profiles
    M = np.array([[0, 0, 0, 0]])
    for x in range(0, trj_len, step):
        print(x)
        S = profile[x]
        F = np.array(S.tolist())
        J = np.insert(F, 0, x, axis=1)
        M = np.append(M, J, axis=0)
    np.savetxt(f"./chain{chain}/cadena{chain}.txt", M)
    H.create_vmd_surface(filename=f"./chain{chain}/hole.vmd", dot_density=10)

def main_f(top,traj,chain_list,chain_ranges,ref_resids,step):
    for chain in chain_list:
        try:
            os.mkdir(os.path.join(f'chain{chain}'))
        except FileExistsError:
            pass

    for chain_resids, chain in zip(chain_ranges, chain_list):
        u = mda.Universe(top, traj)
        chain_traj_xtract("resid " + chain_resids, chain, u)
        run_hole_traj(chain, ref_resids, step)


if __name__ == "__main__":
    CLI = argparse.ArgumentParser()
    
    CLI.add_argument(
        '-p',
        type=str,
        dest='top',
        metavar='Topology',
        help='Topology filename.'
    )

    CLI.add_argument(
        '-x',
        type=str,
        dest='traj',
        metavar='Trajectory',
        help='Trajectory filename.'
    )

    CLI.add_argument(
        '--chain-list',
        nargs='*',
        type=str,
        dest='chain_list',
        metavar='A B C D',
        help= 'List of space separated letters, one for each chain to be analyzed. Example: --chain-list A B'
    )

    CLI.add_argument(
        '--chain-ranges',
        nargs='*',
        type=str,
        dest='chain_ranges',
        metavar='Chain ranges',
        help='List of space separated residue span of each chain to be analyzed. Example: --chain-ranges 1-235 236-470'
    )

    CLI.add_argument(
        '--ref-resids',
        nargs='*',
        type=int,
        dest='ref_resids',
        metavar='Reference residues',
        help='List of space separated numbers, one for each reference residue. Same residues will be used for every chain'
    )

    CLI.add_argument(
        '--step',
        type=int,
        dest='step',
        metavar='Step',
        help='Frame step of hole analysis'
    )

    args = CLI.parse_args()

    main_f(args.top,args.traj,args.chain_list,args.chain_ranges,args.ref_resids,args.step)


    



"""    
    chain_traj_xtract("resid "+chain_range[chain_num],chain_name[chain_num])
    
    
    
for chain_num in range(0,chain_count):
    try:
        os.remove("*old")
        os.remove("*sph")
    except OSError:
        print('Primer Paso')
    hole=run_hole_traj(chain_name[chain_num],ref_resid)
    hole.create_vmd_surface(filename='hole.vmd', dot_density=10)"""
