from ..assistant_system import AssistantSystem
from .funcs_for_importfiles import import_bond_

import os
import numpy as np
import pathlib as Path
import pandas as pd
from tqdm import tqdm
from ase.build import molecule
import time
from collections import deque

import multiprocessing
##################################################
##################################################
class ImportAtoms:
    def __init__(self) -> None:
        self.colm_cell:list = ("cellmin", "cellmax")
        self.colm_mass:list = ("type", "mass")
        self.colm_atom:list = ("id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q")
        
        self.keys_int:list = ("id", "type", "mask")
        
        self.path_wd: Path
        self.systems:tuple[AssistantSystem]
        
        self.cell
        self.mass
        self.atom
        pass
##################################################
    # def import_dumps(self) -> dict:

        
    #     path_input_init:Path.Path = path_inputs[0]
        
    #     with open(path_input_init, "r") as f:
    #         lines = [line.split() for line in f.readlines(1000)]
        
    #     for row, line in enumerate(lines):
    #         if not line:
    #             pass
    #         elif "NUMBER OF ATOMS" in " ".join(line):
    #             atoms_num:int = int(lines[row+1][0])
    #         elif "BOX BOUNDS" in " ".join(line):
    #             cells_num = 3
    #             startrow_cell = row+1
    #         elif "ATOMS id type" in " ".join(line):
    #             startrow_atom = row+1
    #             colm_atom = line[2:]
    #             break
        
    #     self.systems:list = [AssistantSystem() for _ in steps]
        
    #     for 
    #     # cellsize
    #     self.cell:dict = dict()
    #     dates_cell:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_cell, max_rows=cells_num, unpack=True)
    #     for key, date in zip(self.colm_cell, dates_cell):
    #         self.cell[key] = date
    #     # atom
    #     self.atom:dict = dict()
    #     dates_atom:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_atom, max_rows=atoms_num, unpack=True)
    #     sorted_indices = np.argsort(dates_atom[0])
    #     dates_atom = dates_atom[:, sorted_indices]
    #     for key, date in zip(colm_atom, dates_atom):
    #         self.atom[key] = date
    #     if self.colm_atom[2] not in colm_atom:
    #         masks = np.array([0*atoms_num])
    #         self.atom[self.colm_atom[2]] = masks
            
    #     # correction
    #     self.correct_dates()
    #     self.correct_dtype()
    #     # mass
    #     self.make_masses()
    #     return 0
##################################################
##################################################
class ImportBonds:
    def __init__(self) -> None:
        pass
##################################################
    def import_bonds(self):
        path_bond_init = self.path_bonds[0]
        with open(path_bond_init, "r") as f:
            current_row = 0
            while current_row < 5:
                line = f.readline()
                if len(line) == 1:
                    continue
                if 'NUMBER' in line:
                    num_atoms = int(f.readline())
                    current_row += 1 # readline()を行った為
                    skiprows = current_row+1
                    break
                current_row += 1
        
        # 実行プロセスを管理するリスト
        processes = [None]*len(self.path_bonds)
        pool = multiprocessing.Pool(processes=self.num_cores)
        
        for i, path_bond in enumerate(self.path_bonds):            
            # apply_asyncは非同期に関数をpool内で実行する。
            # 関数内で定義した関数を呼び出すことはできないので注意
            processes[i] = pool.apply_async(func=import_bond_, args=(path_bond, skiprows, num_atoms))
            del process
        pool.close()
        
        for system, process in tqdm(zip(self.systems, processes), total=len(processes)):
            system.bondids_tuple, system.bondorders_tuple = process.get()
        pool.join()
        return 0
##################################################
##################################################
##################################################
class ImportFiles(ImportAtoms, ImportBonds):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################