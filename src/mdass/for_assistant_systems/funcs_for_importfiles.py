import pandas as pd
import numpy as np
import itertools
from collections import deque

def import_bond_(path_bond, skiprows, num_atom) -> tuple:
    # 律速
    df_bond:pd.DataFrame = pd.read_csv(path_bond, skiprows=skiprows,  names=[0,1,2,3,4],
                                    sep=" ", header=None, engine="c")
    df_bond:pd.DataFrame = df_bond.replace('Atom', -1).astype("float32")
    array_bond:np.array = df_bond.values
    del df_bond
    
    # bondorderを抽出
    bondids_undent = array_bond[:, 0].astype(int)
    bondorders_undent:np.array = array_bond[:, 1:5]
    
    # 欠損値を含む行を選択するマスクを作成する,bool型をvalueに持つnp.arrayが作成される
    mask_nan = np.any(np.isnan(array_bond), axis=1)
    rownums_mainhead = np.where(mask_nan)[0]
    
    # 欠損値を含む行の2番目と3番目の列を選択し、整数型に変換する
    mainid_bondnum = array_bond[mask_nan, 1:3].astype(int).T
    mainids = mainid_bondnum[0]
    bondnums = mainid_bondnum[1]
    
    # 第二律速
    # 以下2つのarrayは一対一対応でなければならない
    bondorders_tuple = [None]*num_atom
    bondids_tuple = [None]*num_atom
    
    for rownum, mainid, bondnum in zip(rownums_mainhead, mainids, bondnums):
        mainidx = mainid-1 # 以降, main_idxの定義
        bondids_tuple[mainidx] = bondids_undent[rownum+1: rownum+bondnum+1]
        bondorders_tuple[mainidx] = bondorders_undent[rownum+1: rownum+bondnum+1]
    
    bondids_tuple = tuple(bondids_tuple)
    bondorders_tuple = tuple(bondorders_tuple)
    return bondids_tuple, bondorders_tuple
##################################################
def make_connect_tuple_(num_atom:int, bondids_tuple:tuple[int], bondorders_tuple:tuple[list[float]], 
                       cutoff:float=0.5, allow_duplicate:bool=False):
    connect_tuple = [[] for _ in range(num_atom)]
    for mainidx, (bondids_tuple, bondorders_tuple) in enumerate(zip(bondids_tuple, bondorders_tuple)):
        mainid = mainidx+1
        for bondid, bondorder in zip(bondids_tuple, bondorders_tuple):
            if (not allow_duplicate) and (bondid<=mainid):
                continue
            if bondorder[-1] >= cutoff:
                connect_tuple[mainidx].append(bondid)
    connect_tuple:tuple = tuple(connect_tuple)
    return connect_tuple
##################################################
def count_mols_(num_atom:int, connect_tuple:tuple[float], 
                types:np.array, type_elm:dict, 
                min_mol:int=1, max_mol:int=100000, condition=None):
    if condition is None:
        target_atoms = np.array([True]*num_atom)
    else:
        target_atoms = target_atoms[condition]
    
    type_set = set(types)
    mols_temp:dict = dict()
    
    #幅優先探索
    visited_ids = [0]*num_atom
    for mainidx in range(num_atom):
        if visited_ids[mainidx]:
            continue
        if not target_atoms[mainidx]:
            continue
        current_mol:list = [0]*max(type_set)
        
        que = deque([mainidx])
        while que:
            current_idx = que.popleft()
            if visited_ids[current_idx]:
                continue
            visited_ids[current_idx] = 1
            current_mol[types[current_idx]-1] += 1
            for add_id in connect_tuple[current_idx]:
                add_idx = add_id-1
                if visited_ids[add_idx]:
                    continue
                que.append(add_idx)
        
        if min_mol <= sum(current_mol) <= max_mol:
            current_mol_tuple:tuple = tuple(current_mol)
            if current_mol_tuple not in mols_temp:
                mols_temp[current_mol_tuple] = 0
            mols_temp[current_mol_tuple] += 1
    mols = dict()
    for mol_type, amount in mols_temp.items():
        mol_elm:list = [f"{type_elm[i+1]}{stoichiometry} " for i,stoichiometry in enumerate(mol_type) if stoichiometry!=0]
        mol_str:str = "".join(mol_elm)
        mols[mol_str] = amount
    return mols
##################################################
def count_trios_(num_atom:int, connect_tuple:tuple[int], types:np.array, 
                 combis_type, type_elm, condition=None):
    if condition is None:
        target_atoms = np.array([True]*num_atom)
    else:
        target_atoms = target_atoms[condition]
    
    trios_temp:dict = {combi:0 for combi in combis_type}
    for mainidx, bondids in enumerate(connect_tuple):
        if not target_atoms[mainidx]:
            continue
        maintype = types[mainidx]
        for bondid1, bondid2 in itertools.combinations(bondids, 2):
            if (bondid1<=mainidx+1) or (bondid2<=mainidx+1):
                continue
            bondtype1 = types[bondid1-1]
            bondtype2 = types[bondid2-1]
            sorted_list = sorted([bondtype1, bondtype2])
            sorted_list.insert(1, maintype)
            combi_now = tuple(sorted_list)
            trios_temp[combi_now] += 1
    
    trios:dict = dict()
    for combi_type, amount in trios_temp.items():
        combi_elm:list = [type_elm[type_] for type_ in combi_type]
        combi_str:str = "-".join(combi_elm)
        trios[combi_str] = amount
    return trios