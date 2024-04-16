from .assistant_system import AssistantSystem
from .for_assistant_systems.import_files import ImportFiles
from .for_assistant_systems.analyses import Analyses

# from .dataset import mass_dict

import os
import pathlib as Path
import multiprocessing
##################################################
##################################################
class AssistantSystems(AssistantSystem,
                       ImportFiles,
                       Analyses,
                       ):
    def __init__(self, 
                 path_wd:"str"="./",
                 ):
        self.path_wd:Path.Path = Path.Path(path_wd)
        self.systems:tuple[AssistantSystem] = tuple()
        
        self.colm_cell:list = ("cellmin", "cellmax")
        self.colm_mass:list = ("type", "mass")
        self.colm_atom:list = ("id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q")
        
        self.keys_int:list = ("id", "type", "mask")
        
        self.num_cores = multiprocessing.cpu_count()
        print("Number of CPU cores:",self.num_cores)
##################################################
    def get_filename_pos(self):
        path_wd:Path.Path = self.path_wd
        steps_pos:list[str] = []
        for name_file in os.listdir(path_wd):
            if "dump.pos." in name_file:
                steps_pos.append(int(name_file[9:]))
        steps_pos.sort()
        self.path_poses = tuple(self.path_wd.joinpath(f"dump.pos.{step}") for step in steps_pos)
        
        pos_filemum = len(self.path_poses)
        if  self.systems:
            assert len(self.systems)==pos_filemum, "読み込むファイル数は統一する必要がある"
        else:
            self.systems = tuple(AssistantSystem(self.path_wd) for _ in range(pos_filemum))
##################################################
    def get_filename_bond(self):
        path_wd:Path.Path = self.path_wd
        self.steps_bond:list[str] = []
        for name_file in os.listdir(path_wd):
            if "dump.bond." in name_file:
                self.steps_bond.append(int(name_file[10:]))
        self.steps_bond.sort()
        self.path_bonds = tuple(self.path_wd.joinpath(f"dump.bond.{step}") for step in self.steps_bond)
        
        bond_filemum = len(self.path_bonds)
        if  self.systems:
            assert len(self.systems)==bond_filemum, "読み込むファイル数は統一する必要がある"
        else:
            self.systems = tuple(AssistantSystem(self.path_wd) for _ in range(bond_filemum))
##################################################
##################################################
##################################################