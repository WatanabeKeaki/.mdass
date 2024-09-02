import numpy as np
from collections import deque
import itertools
import multiprocessing
from tqdm import tqdm

from .funcs_for_importfiles import make_connect_tuple_
from .funcs_for_importfiles import count_mols_
from .funcs_for_importfiles import count_trios_
##################################################
##################################################
class AnalyseBonds:
    def __init__(self) -> None:
        pass
##################################################
    def make_connect_tuple_forsystems(self, cutoff:float=0.5, allow_duplicate:bool = False):
        num_atom = self.atom["id"].shape[0]
        processes = [None]*len(self.path_bonds)
        pool = multiprocessing.Pool(processes=self.num_cores)
        
        for i, system in enumerate(self.systems):
            processes[i] = pool.apply_async(func=make_connect_tuple_, 
                                            args=(num_atom, 
                                                  system.bondids_tuple, 
                                                  system.bondorders_tuple, 
                                                  cutoff, 
                                                  allow_duplicate)
                                            )
        pool.close()
        for system, process in tqdm(zip(self.systems, processes), total=len(processes)):
            system.connect_tuple = process.get()
        pool.join()
        return 0
##################################################
    def count_mols_forsystems(self, min_mol:int=1, max_mol:int=100000, condition=None):
        
        
        num_atom:int = self.atom["id"].shape[0]
        processes = [None]*len(self.path_bonds)
        pool = multiprocessing.Pool(processes=self.num_cores)
        types = self.atom["type"]
        type_elm = self.type_elm
        
        for i, system in enumerate(self.systems):
            processes[i] = pool.apply_async(func=count_mols_, 
                                            args=(num_atom, 
                                                  system.connect_tuple, 
                                                  types, 
                                                  type_elm, 
                                                  min_mol, 
                                                  max_mol, 
                                                  condition)
                                            )
        pool.close()
        
        molsdicts_temp = {step: None for step in self.steps_bond}
        for step, process in tqdm(zip(self.steps_bond, processes), total=len(processes)):
            molsdicts_temp[step] = process.get()
        pool.join()
        
        mol_set:set = set()
        for mols_dict in molsdicts_temp.values():
            for mol in mols_dict.keys():
                mol_set.add(mol)
        
        amounts = [0]*len(self.steps_bond)
        self.mols_dicts = {mol: amounts.copy() for mol in mol_set}
        
        for i, molsdict_temp in enumerate(molsdicts_temp.values()):
            for mol, amount in molsdict_temp.items():
                self.mols_dicts[mol][i] = amount
        return 0
##################################################
    def count_trios_forsystems(self, condition=None) -> int:
        bondlength = 3
        
        num_atom:int = self.atom["id"].shape[0]
        processes = [None]*len(self.path_bonds)
        pool = multiprocessing.Pool(processes=self.num_cores)
        types = self.atom["type"]
        combis_type = self.generate_unique_combis(bondlength)
        type_elm = self.type_elm
        
        for i, system in enumerate(self.systems):
            processes[i] = pool.apply_async(func=count_trios_, 
                                            args=(num_atom, 
                                                  system.connect_tuple, 
                                                  types, 
                                                  combis_type, 
                                                  type_elm, 
                                                  condition)
                                            )
        pool.close()
        triosdicts_temp0 = {step: None for step in self.steps_bond}
        for step, process in tqdm(zip(self.steps_bond, processes), total=len(processes)):
            triosdicts_temp0[step] = process.get()
        pool.join()
        
        combi_set:list = [combi for combi in triosdicts_temp0[self.steps_bond[0]].keys()]
        amounts = [0]*len(self.steps_bond)
        trios_dicts_temp1 = {combi: amounts.copy() for combi in combi_set}
        for i, triosdict_temp in enumerate(triosdicts_temp0.values()):
            for combi, amount in triosdict_temp.items():
                trios_dicts_temp1[combi][i] = amount
        
        self.trios_dicts = dict()
        for combi, amounts in trios_dicts_temp1.items():
            if not all(amount==0 for amount in amounts):
                self.trios_dicts[combi] = amounts
        
        return 0
##################################################
    def count_duet(self):
        pass
##################################################
    def count_quartet(self):
        pass
##################################################
##################################################
##################################################
class Analyses(AnalyseBonds):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################