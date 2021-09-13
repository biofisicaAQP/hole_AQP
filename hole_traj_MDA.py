import numpy as np
import MDAnalysis as mda
import MDAnalysis.analysis.hole
import MDAnalysis.transformations.translate
from MDAnalysis.analysis.hole import HOLEtraj
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

def chain_traj_xtract(resid_range,chain):
    prot=u.select_atoms(resid_range)
    prot.write("{}coord.pdb".format(chain))
    prot.write("{}traj.xtc".format(chain),frames="all")
#    with MDAnalysis.Writer("{}traj.xtc".format(chain), prot.n_atoms) as W:
#        for ts in u.trajectory:
#            W.write(prot)     

def chain_rms_fit(chain,ref_res):
    raw_traj = "{}traj.xtc".format(chain)
    coord= "{}coord.pdb".format(chain)
    ref_frame = mda.Universe(coord,coord)
    raw_univ = mda.Universe(coord,raw_traj)
    alignment = align.AlignTraj(raw_univ, ref_frame, filename="{}rmsfit.dcd".format(chain))
    alignment.run()
    align_traj = "{}rmsfit.dcd".format(chain)
    align_univ = mda.Universe(coord,align_traj)
    refAtoms = align_univ.residues[ref_res-u.residues[0].resid].atoms
    center_transformation = MDAnalysis.transformations.center_in_box(refAtoms, center='mass', point=[0,0,0])
    align_univ.trajectory.add_transformations(center_transformation)
    return(align_univ)

def run_hole_traj(chain,ref_res):
    fit=chain_rms_fit(chain,ref_res)              
    H = HOLEtraj(fit, executable="~/hole2/exe/hole",endrad=5.0,cpoint=True,cvect=[0, 0, -1], sample=0.125)
    H.run()
    profile = H.profiles
    M = np.array([[0,0,0,0]])
    for  x in range(0,len(fit.trajectory)):
        S = profile[x]
        F = np.array(S.tolist())
        J = np.insert(F, 0, x, axis=1)
        M = np.append(M, J, axis=0)
    np.savetxt("cadena{}.txt".format(chain), M)

COR = "top.pdb" #sustituir por el nombre de la topología de la proteína (formato PDB)
TRJ = "../pbcprot_30ns.dcd" #sustituir por el nombre de la trayectoria (formato DCD)
u = mda.Universe(COR,TRJ)
chain_name=["A","B","C","D"] #definir el nombre de las cadenas (si se va a medir poro central definir 1 solo elemento en esta lista)
chain_range=["2-245","247-490","492-737","739-982"] #definir el rango de residuos para cada cadena (si se va a medir poro central, definir el rango total de residuos)
ref_resid= 200 #definir el residuo que se va a usar como origen de coordenadas. Conviene usar un residuo que sea más o menos rígido respecto de la estructura.
chain_count=len(chain_name)


for chain_num in range(0,chain_count):
    chain_traj_xtract("resid "+chain_range[chain_num],chain_name[chain_num])

for chain_num in range(0,chain_count):    
    run_hole_traj(chain_name[chain_num],ref_resid)
