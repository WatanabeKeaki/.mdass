import numpy as np
from collections import deque
import itertools
##################################################
##################################################
class AnalyseBond:
    def __init__(self) -> None:
        pass
##################################################
    def make_connect_list(self, cutoff:float=0.5, allow_duplicate:bool = False):
        num_atom = self.atom["id"].shape[0]
        self.connect_list = [[] for _ in range(num_atom)]
        for mainidx, (bondids, bondorders) in enumerate(zip(self.bondids_list, self.bondorders_list)):
            mainid = mainidx+1
            for bondid, bondorder in zip(bondids, bondorders):
                if (not allow_duplicate) and (bondid<=mainid):
                    continue
                if bondorder[-1] >= cutoff:
                    self.connect_list[mainidx].append(bondid)
##################################################
    def count_mols(self, min_mol:int=1, max_mol:int=100000, condition=None):
        num_atom = self.atom["id"].shape[0]
        self.connect_list:list
        
        # 探索する原子の条件、type==3など。
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        types = self.atom["type"]
        type_set = set(types)
        
        self.mols:dict = dict()
        
        #幅優先探索
        visited_ids = [0]*num_atom
        for mainidx in range(num_atom):
            if visited_ids[mainidx]:
                continue
            if not target_atoms[mainidx]:
                continue
            current_mol = [0]*max(type_set)
            
            que = deque([mainidx])
            while que:
                current_idx = que.popleft()
                if visited_ids[current_idx]:
                    continue
                visited_ids[current_idx] = 1
                current_mol[types[current_idx]-1] += 1
                for add_id in self.connect_list[current_idx]:
                    add_idx = add_id-1
                    if visited_ids[add_idx]:
                        continue
                    que.append(add_idx)
            
            if min_mol <= sum(current_mol) <= max_mol:
                current_mol_tuple = tuple(current_mol)
                if current_mol_tuple not in self.mols:
                    self.mols[current_mol_tuple] = 0
                self.mols[current_mol_tuple] += 1
##################################################
    def count_trios(self, condition=None):
        bondlength:int = 3
        
        num_atom = self.atom["id"].shape[0]
        self.connect_list:list
        
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        types = self.atom["type"]
        combis = self.generate_unique_combis(bondlength)
        
        trios_temp:dict = {combi:0 for combi in combis}
        for mainidx, bondids in enumerate(self.connect_list):
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
        self.trios = {key: value for key, value in trios_temp.items()}
##################################################
    def count_duet(self):
        pass
##################################################
    def count_quartet(self):
        pass
##################################################
##################################################
##################################################
class Analyse(AnalyseBond):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################