import MDAnalysis.transformations.translate, os, numpy as np, MDAnalysis as mda
from MDAnalysis.analysis import hole2
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd


def chain_traj_xtract(resid_range,chain):
    prot=u.select_atoms(resid_range)
    prot.write("{}coord.pdb".format(chain))
    with MDAnalysis.Writer("{}traj.xtc".format(chain), prot.n_atoms) as W:
        for ts in u.trajectory:
            W.write(prot)     

def chain_rms_fit(chain,ref_res):
    raw_traj = "{}traj.xtc".format(chain)
    coord= "{}coord.pdb".format(chain)
    ref_frame = mda.Universe(coord,coord)
    raw_univ = mda.Universe(coord,raw_traj)
    alignment = align.AlignTraj(raw_univ, ref_frame, filename="{}rmsfit.dcd".format(chain))
    alignment.run()
    align_traj = "{}rmsfit.dcd".format(chain)
    align_univ = mda.Universe(coord,align_traj)
    refAtoms = align_univ.residues[ref_res-1].atoms
    center_transformation = MDAnalysis.transformations.center_in_box(refAtoms, center='mass', point=[0,0,0])
    align_univ.trajectory.add_transformations(center_transformation)
    return(align_univ)

def run_hole_traj(chain,ref_res):
    fit=chain_rms_fit(chain,ref_res)              
    H = hole2.HoleAnalysis(fit, executable="~/hole2/exe/hole")
    H.run()
    profile = H.profiles
    M = np.array([[0,0,0,0]])
    for  x in range(0,len(fit.trajectory)):
        S = profile[x]
        F = np.array(S.tolist())
        J = np.insert(F, 0, x, axis=1)
        M = np.append(M, J, axis=0)
    np.savetxt("cadena{}.txt".format(chain), M)

COR = "first.pdb" #sustituir por el nombre de la topología de la proteína (formato PDB)
TRJ = "top.dcd" #sustituir por el nombre de la trayectoria (formato DCD)
u = mda.Universe(COR,TRJ)
chain_name=["X"] #definir el nombre de las cadenas (si se va a medir poro central definir 1 solo elemento en esta lista)
chain_range=["1-1140"] #definir el rango de residuos para cada cadena (si se va a medir poro central, definir el rango total de residuos)
ref_resid= 140 #definir el residuo que se va a usar como origen de coordenadas. Conviene usar un residuo que sea más o menos rígido respecto de la estructura.
chain_count=len(chain_name)


for chain_num in range(0,chain_count):
    chain_traj_xtract("resid "+chain_range[chain_num],chain_name[chain_num])

for chain_num in range(0,chain_count):
    try:
        os.remove("*old*")
    except OSError:
        print('Primer_Paso')
    run_hole_traj(chain_name[chain_num],ref_resid)


