import pf_calculator.pf_calculator as pf_calc
import os
import matplotlib.pyplot as plt
import pandas as pd

molecule = os.path.dirname(os.getcwd()).split('/')[-1]
dirs = next(os.walk('.'))[1]
dirs.sort()
chains = [0,1,2,3]
chains_letter = ['A','B','C','D']
timestep = 10**(-12)
(ref_z_1, ref_z_2, ref_xy_1, ref_xy_2) = (0,1,1,2)
pf_dict_list = []
for dir in dirs:
    filename = os.path.join(dir, 'perm.nc')
    for chain_id in chains:
        (pf,plot,compendio_atomos) = pf_calc.pf_calculator(filename, chain_id, ref_z_1, ref_z_2, ref_xy_1, ref_xy_2, timestep=timestep)
        plt.savefig(f'{dir}_{chains_letter[chain_id]}_msd.png')
        plt.close()
        print(f'{dir}_{chains_letter[chain_id]}_pf: {pf/1e-14:.4f}e-14')
        pf_dict = {'Time':dir,'Molecule': molecule ,'Chain':chains_letter[chain_id], 'pf':pf}
        pf_dict_list.append(pf_dict)
        fname_atom_list =f'atoms_in_pore{chains_letter[chain_id]}_{dir}.txt'
        with open(fname_atom_list,'w') as f:
            for atom in compendio_atomos:
                f.write("%s\n" % atom)

pf_dataframe = pd.DataFrame(pf_dict_list)
print(pf_dataframe)
pf_dataframe.to_csv('pf_dataframe.csv')

