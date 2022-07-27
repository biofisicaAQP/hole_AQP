#!/usr/bin/python
import argparse
from turtle import st
import matplotlib.pyplot as plt
from pmf_calculator.pmf_calculator import pmf_calculator as pmf_calc
from pmf_calculator.pmf_calculator import histogram_calculator as histogram_calc
from lib.functions import *

parser = argparse.ArgumentParser()
parser.add_argument(
    '-x',
    nargs='*',
    type=str,
    dest= 'trajs',
    metavar='TRAJS list',
    help='List of space separated trajectory filenames. Example: -x wat.nc aox.nc'
)

parser.add_argument(
    '-m',
    nargs='*',
    type=str,
    dest= 'molecules',
    metavar='MOL',
    help='List of space separated molecule names. Example: -m WAT AOX'
)

parser.add_argument(
    '--chain-list',
    nargs='*',
    type=str,
    dest='chain_list',
    metavar='CHAIN',
    help= 'List of space separated letters, one for each chain to be analyzed. Example: --chain-list A B'
)

parser.add_argument(
    '-z',
    nargs=2,
    type=int,
    dest= 'z_ref',
    metavar='zREF',
    help='List of space separated z reference atoms. Example: -z 0 1'
)

parser.add_argument(
    '-xy',
    nargs=2,
    type=int,
    dest= 'xy_ref',
    metavar='xyREF',
    help='List of space separated xy reference atoms. Example: -xy 1 2'
)

args = parser.parse_args()
file_list = args.trajs
molecule_list = args.molecules
chains_letter = args.chain_list
(ref_z_1, ref_z_2) = args.z_ref
(ref_xy_1, ref_xy_2) = args.xy_ref

histogram_dict_list =[]
for chain_id, chain in enumerate(chains_letter):
    for filename, molecule in zip(file_list, molecule_list, ):
        mol_dict = {}
        mol_dict['histogram'],mol_dict['hist_edges'],mol_dict['array_shape'] = histogram_calc(filename,chain_id, ref_z_1, ref_z_2, ref_xy_1, ref_xy_2,bins=1000, pore_radius=4)
        mol_dict['molecule'] = molecule
        mol_dict['chain'] = chain
        histogram_dict_list.append(mol_dict)
total_histogram_dict = {}

for chain in chains_letter:
    total_histogram = 0
    for entry in histogram_dict_list:
        if entry['chain'] == chain:
            total_histogram += entry['histogram']
    total_histogram_dict[chain] = total_histogram
    for entry in histogram_dict_list:
        if entry['chain'] == chain:
            entry['relative_histogram'] = entry['histogram']/total_histogram

for entry in histogram_dict_list:
    entry['pmf'], entry['edges_pmf'] = pmf_calc(entry['histogram'],entry['hist_edges'], entry['array_shape'])
    chain = entry['chain']
    molecule = entry['molecule']
    plt.plot(entry['edges_pmf'],entry['pmf'])
    plt.title(chain + ' ' + molecule)
    plt.xlim(20,120)
    plt.savefig(f'chain_{chain}_{molecule}_pmf.png')
    plt.clf()

for chain in chains_letter:
    fig,axs = plt.subplots(len(molecule_list),1)
    fig.suptitle(f'Chain {chain}')
    mol_counter = 0
    for entry in histogram_dict_list:
        if chain == entry['chain']:
            start = len(entry['edges_pmf'])//10
            stop = len(entry['edges_pmf']) - start
            axs[mol_counter].plot(entry['edges_pmf'][start:stop],entry['pmf'][start:stop])
            axs[mol_counter].set_title(entry['molecule'])
#            axs[mol_counter].set_ylim(0,10)
            mol_counter+=1
    plt.setp(axs, ylim=(2,10))
    plt.savefig(f'chain_{chain}_pmf.png')
    plt.clf()