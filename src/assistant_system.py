import import_file
import correction
import pathlib as Path
from tqdm import tqdm
import pandas as pd
import dateset
from time import time

class AssistantSystem(import_file.ImportAtom,
                      import_file.ImportPara,
                      import_file.ImportBond,   
                      correction.CorrectAtom):
    def __init__(self, path_wd:"str"="./"):
        self.path_wd:Path.Path = Path.Path(path_wd)
        
        self.colm_cell:list = ["cellmin", "cellmax"]
        self.colm_mass:list = ["type", "mass"]
        self.colm_atom:list = ["id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q"]
        
        self.keys_int:list = ["id", "type", "mask"]
        
        self.cell = None
        self.masses = None
        self.atoms = None
        pass

test = AssistantSystem()
# test.import_rd_test("../input.rd")
test.import_para("../../testfiles/para.rd")
# test.import_bond("../test/dump.bond.00")

# time0 = time()
for _ in tqdm(range(100)):
    # test.import_rd("../../testfiles/input.rd")
    test.import_bond("../../testfiles/dump.bond.0")
    # test.import_dump("../../testfiles/dump.pos.0")
    # test.import_para("../../testfiles/para.rd")
# time1 = time()
# print(time1-time0)
# test.import_dump("../dump.pos.0")