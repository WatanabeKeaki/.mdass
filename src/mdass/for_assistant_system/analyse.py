import numpy as np
from collections import deque
import itertools
# from .cppfiles.make_connect_list import make_connect_list
##################################################
##################################################
class AnalyseBond:
    def __init__(self) -> None:
        pass
##################################################
    def make_connectids_list(self, cutoff:float=0.5, allow_duplicate:bool = True):
        """
        allow_duplicateまわりにバグがありそうなので、
        Falseを使わないこと
        """
        num_atom = self.atom["id"].shape[0]
        self.connectids_list:list = [[] for _ in range(num_atom)]
        for mainidx, (bondids, bondorders) in enumerate(zip(self.bondids_list, self.bondorders_list)):
            mainid = mainidx+1
            for bondid, bondorder in zip(bondids, bondorders):
                # if (not allow_duplicate) and (bondid<=mainid):
                #     continue
                if bondorder[-1] >= cutoff:
                    self.connectids_list[mainidx].append(bondid)
        # self.connectids_list:list = make_connect_list(num_atom, 
        #                                               self.bondids_list, self.bondorders_list, 
        #                                               cutoff, allow_duplicate)
        return 0
##################################################
    def count_mols(self, min_mol:int=1, max_mol:int=100000, condition=None):
        num_atom = self.atom["id"].shape[0]
        self.connectids_list:list
        
        # 探索する原子の条件、type==3など。
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        types = self.atom["type"]
        type_set = set(types)
        
        self.mols:dict = dict()
        
        #幅優先探索
        visited_idxs = [0]*num_atom
        for mainidx in range(num_atom):
            if visited_idxs[mainidx]:
                continue
            if not target_atoms[mainidx]:
                continue
            
            current_mol = [0]*max(type_set)
            
            que = deque([mainidx])
            while que:
                current_idx = que.popleft()
                if visited_idxs[current_idx]:
                    continue
                visited_idxs[current_idx] = 1
                current_mol[types[current_idx]-1] += 1
                for add_id in self.connectids_list[current_idx]:
                    add_idx = add_id-1
                    if visited_idxs[add_idx]:
                        continue
                    que.append(add_idx)
            
            if min_mol <= sum(current_mol) <= max_mol:
                current_mol_tuple = tuple(current_mol)
                if current_mol_tuple not in self.mols:
                    self.mols[current_mol_tuple] = 0
                self.mols[current_mol_tuple] += 1
        return 0
##################################################
    def count_pair(self, condition=None):
        bondlength:int = 2
        
        num_atom = self.atom["id"].shape[0]
        self.connectids_list:list
        
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        types = self.atom["type"]
        combis = self.generate_unique_combis(bondlength)
        
        pairs_temp:dict = {combi:0 for combi in combis}
        for mainidx, bondids in enumerate(self.connectids_list):
            if not target_atoms[mainidx]:
                continue
            maintype = types[mainidx]
            for bondid in bondids:
                if (bondid<=mainidx+1):
                    continue
                bondtype = types[bondid-1]
                sorted_list = sorted([maintype, bondtype])
                combi_now = tuple(sorted_list)
                pairs_temp[combi_now] += 1
        self.pairs = {key: value for key, value in pairs_temp.items()}
        return 0
##################################################
    def count_trios(self, condition=None):
        bondlength:int = 3
        
        num_atom = self.atom["id"].shape[0]
        self.connectids_list:list
        
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        types = self.atom["type"]
        combis = self.generate_unique_combis(bondlength)
        
        trios_temp:dict = {combi:0 for combi in combis}
        for mainidx, bondids in enumerate(self.connectids_list):
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
        return 0
##################################################
    def count_duet(self):
        pass
##################################################
    def count_quartet(self):
        pass
##################################################
    def delete_mols(self, min_mol:int=5, condition=None):
        num_atom = self.atom["id"].shape[0]
        self.connectids_list:list
        
        # 探索する原子の条件、type==3など。
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        #幅優先探索
        visited_idxs:list = [0]*num_atom
        remain_idxs:list = [0]*num_atom
        for mainidx in range(num_atom):
            if visited_idxs[mainidx]:
                continue
            if not target_atoms[mainidx]:
                continue
            
            current_mol_idxs:list = []
            
            que = deque([mainidx])
            while que:
                current_idx = que.popleft()
                if visited_idxs[current_idx]:
                    continue
                visited_idxs[current_idx] = 1
                current_mol_idxs.append(current_idx)
                for add_id in self.connectids_list[current_idx]:
                    add_idx = add_id-1
                    if visited_idxs[add_idx]:
                        continue
                    que.append(add_idx)
                
            if min_mol <= len(current_mol_idxs):
                # print(current_mol_idxs)
                for idx in current_mol_idxs:
                    remain_idxs[idx] = 1
        
        # print(remain_idxs)
        # exit()
        # あとでこれ以下の部分は別の関数にして分離
        temp_atom = dict()
        remain_idxs = np.array(remain_idxs)
        for key, value in self.atom.items():
            temp_atom[key] = value.copy()[remain_idxs==1]
        self.atom = temp_atom.copy()
        return 0
##################################################
    def delete_bulks(self, max_mol:int=100, condition=None):
        num_atom = self.atom["id"].shape[0]
        print(num_atom)
        self.connectids_list:list
        
        # 探索する原子の条件、type==3など。
        if condition is None:
            target_atoms = np.array([True]*num_atom)
        else:
            target_atoms = target_atoms[condition]
        
        #幅優先探索
        visited_idxs:list = [0]*num_atom
        remain_idxs:list = [0]*num_atom
        for mainidx in range(num_atom):
            if visited_idxs[mainidx]:
                continue
            if not target_atoms[mainidx]:
                continue
            
            current_mol_idxs:list = []
            
            que = deque([mainidx])
            while que:
                current_idx = que.popleft()
                if visited_idxs[current_idx]:
                    continue
                visited_idxs[current_idx] = 1
                current_mol_idxs.append(current_idx)
                for add_id in self.connectids_list[current_idx]:
                    add_idx = add_id-1
                    if visited_idxs[add_idx]:
                        continue
                    que.append(add_idx)
                
            if len(current_mol_idxs) <= max_mol:
                # print(current_mol_idxs)
                for idx in current_mol_idxs:
                    remain_idxs[idx] = 1
        
        # print(remain_idxs)
        # exit()
        # あとでこれ以下の部分は別の関数にして分離
        temp_atom = dict()
        remain_idxs = np.array(remain_idxs)
        for key, value in self.atom.items():
            temp_atom[key] = value.copy()[remain_idxs==1]
        self.atom = temp_atom.copy()
        return 0
##################################################
    def make_neighbor_list(self):
        """
        1. セルをメッシュで分ける
        2. メッシュに属する原子idを格納したリストを作成
        3. あるメッシュがどのメッシュと隣接しているかを判定
        """
        
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