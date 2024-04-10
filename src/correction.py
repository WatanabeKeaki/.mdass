import dateset
import numpy as np
from ase.build import molecule

class CorrectAtom:
    def __init__(self) -> None:
        self.colm_cell:list = ["cellmin", "cellmax"]
        self.colm_mass:list = ["type", "mass"]
        self.colm_atom:list = ["id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q"]
        
        self.keys_int:list = ["id", "type", "mask"]
        pass
##################################################
    def correct_dates(self):
        """_summary_
        cellsizeの最小値が負の場合、laichはerrorを発生する。\n
        それ以外も解析をしやすくする目的でcellsizeの最小値を0に合わせる。\n
        それに応じてcellsizeの最大値や原子のxyz座標も調整する。\n
        inputfileを読み込み、cellsizeやatomsをinstanceが持っていれば実行可能。
        """
        cell_min:np.array = self.cell[self.colm_cell[0]].copy()
        # positions
        for i, colm in enumerate(self.colm_atom[3:6]):
            self.atoms[colm] -= cell_min[i]
        # cellsizes
        for colm in self.colm_cell:
            self.cell[colm] -= cell_min
        del cell_min
        return 0
##################################################
    def correct_dtype(self):
        """_summary_
        atomsのid, type, maskはnp.int32型である必要がある。\n
        仕様上、読み込まれたままのデータは上記3つがnp.float64型なので、\n
        intに直す。
        """
        for colm in self.colm_atom[0:3]:
            self.atoms[colm] = self.atoms[colm].astype(np.int64)
        return 0
##################################################
    def make_masses(self):
        """_summary_
        self.massesを作成する。
        self.atomsとself.elm_typeが存在する事が条件。
        よって、inputやdumoposと、para.rdを読み込んでいる必要がある。
        self.elm_typeを手動で{1,"C", 2:"H"}のように設定も可能。
        """
        mass_dict:dict = dateset.mass_dict # from dateset.py
        types_set = set(self.atoms[self.colm_atom[1]])
        self.masses = {
        self.colm_mass[0]: np.array([type_ for type_ in self.type_elm.keys() if type_ in types_set]), 
        self.colm_mass[1]: np.array([mass_dict[element] for type_, element in self.type_elm.items() if type_ in types_set])
        }
        return 0
##################################################
    def sort_by_id(self):
        """_summary_
        idの順番がバラバラになっている時、atomsをid順に直す。
        """
        sorted_index = np.argsort(self.atoms['id'])
        temp_atoms = dict()
        for key, values in self.atoms.items():
            temp_atoms[key] = [values[i] for i in sorted_index]
        self.atoms = temp_atoms
        del sorted_index
        del temp_atoms
        return 0