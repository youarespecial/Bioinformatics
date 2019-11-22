import os
import sys
import glob
import time
import argparse
import pandas as pd
import numpy as np

def get_ref_pair(ref_path):
    pair_set = set()
    with open(ref_path) as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            sp = line.strip().split()
            if len(sp) != 2:
                print('Warn split error, len is not 2 : ', sp)
                continue
            pair = tuple(sp)
            if not pair in pair_set:
                pair_set.add(pair)
            else:
                print('Warning found same pair: ', pair)
    return pair_set

def get_fname_pair(input_dir):
    result = {}
    fname_li = os.listdir(input_dir)
    left_fnames = [i for i in fname_li if i.endswith('_L.txt')]
    right_fnames = [i for i in fname_li if i.endswith('_R.txt')]
    pair_li = [(left, right) for left in left_fnames for right in right_fnames]
    return pair_li

def read_gene(path):
    gene_li = []
    with open(pair) as f:
        for idx, line in enumerate(f):
            if idx == 0:
                continue
            sp = line.strip().split()
            if len(sp) != 2:
                print('read gene split error: ', sp)
                continue
            gene = sp[1].strip().strip('"')
            gene_li.append(gene)
    return gene_li

def get_paired_gene(fname_pair, input_dir):
    pair_set = set()
    paths = [os.path.join(input_dir, fname ) for fname in fname_pair]
    left_genes, right_genes = [read_gene(p) for p in paths]
    nrof_left, nrof_right = len(left_genes), len(right_genes)
    if nrof_left and nrof_right:
        pair_li = [(left, right) for left in left_genes for right in right_genes]   
        pair_set = set(pair_li)
    return pair_set

def do_gene_query(txt_dir, ref_path, output_dir):
    ref_pairs = get_ref_pair(ref_path)
    nrof_ref_pair = len(ref_pairs)
    print('Found %d pairs in ref path'%(nrof_ref_pair))

    fname_pair_li = get_fname_pair(txt_dir)
    out_li = []
    header = ['Left-Fname', 'Right-Fname', 'Left-Gene', 'Right-Gene']
    nrof_file_pair = len(fname_pair_li)
    for idx, fname_pair in enumerate(fname_pair_li, start=1):
        pair_set = get_paired_gene(fname_pair, txt_dir)
        common = pair_set & nrof_ref_pair
        if len(common):
            for gene_pair in common:
                print('Process {:3d}/{:3d} file pair, found gene pair: {} in file pair: {}'.format(
                    idx, nrof_file_pair, gene_pair, fname_pair))
                li = list(fname_pair) +list(gene_pair)
                out_li.append(li)
    if out_li:
        df = pd.DataFrame(out_li, columns=header, index=None)
        print(df.head(10))
        csv_path = os.path.join(output_dir, 'gene_query.csv')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        df.to_csv(csv_path, index=None, sep='\t')
        print('Saving to ', csv_path)
    else:
        print('Found no paired gene in all %d paired files'%(nrof_file_pair))


    print('All done!!!')

if __name__ == '__main__':
    input_dir = '/home/wzk/Data/hfz/RL'
    ref_path = '/home/wzk/Data/hfz/RL/ReceptorLigand.txt'
    output_dir = './'
    do_gene_query(input_dir, ref_path, output_dir)
