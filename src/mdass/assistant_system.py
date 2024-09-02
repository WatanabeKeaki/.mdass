from .for_assistant_system.import_file import ImportFile
from .for_assistant_system.analyse import Analyse
from .for_assistant_system.export_file import ExportFile
from .for_assistant_system.exec import Exec
from .for_assistant_system.operate import Operate
from .for_assistant_system.select_atom import SelectAtom

from .dataset import mass_dict

import pathlib as Path
import numpy as np
from tqdm import tqdm
from time import time
import itertools
##################################################
##################################################
class AssistantSystem(ImportFile,
                      Analyse,
                      ExportFile, 
                      Exec, 
                      Operate, 
                      SelectAtom, 
                       ):
    def __init__(self, path_wd:"str"="./"):
        self.path_wd:Path.Path = Path.Path(path_wd)
        
        self.colm_cell:tuple = ("min", "max")
        self.colm_mass:tuple = ("type", "mass")
        self.colm_atom:tuple = ("id", "type", "mask", 
                                "x", "y", "z", 
                                "vx", "vy", "vz", 
                                "fx", "fy", "fz", 
                                "q")
        self.colm_press:tuple = ("Step", 
                                 "CellX", "CellY", "CellZ", 
                                 "StressXX", "StressYY", "StressZZ", 
                                 "StressXY", "StressYZ", "StressZX")
        
        self.keys_int:tuple = ("id", "type", "mask")
        
        self.cell:dict = dict()
        self.mass:dict = dict()
        self.atom:dict = dict()
        self.press:dict = dict()
        
        self.keyword_groups = (
                              ('Mode', 'XYZFile', 'ParaFile', 'Unit', 'ForceField'), 
                              ('TimeStep', 'TotalStep'), 
                              ('ObserveStep', 'FileStep', 'SeparateDump', 'BondStep', 'StressStep', 
                               "SaveRestartStep", "DumpBondsOvito", "BondOrderCUTOFF"), 
                              ('AtomPress', 'AtomForce', "ShowMask"),
                              ('MPIGridX', 'MPIGridY', 'MPIGridZ'), 
                              ('OMPGridX', 'OMPGridY', 'OMPGridZ'), 
                              ('DelR', 'MaxR', 'Factor'), 
                              ('CUTOFF', 'MARGIN', 'GhostFactor'), 
                              ('Conv_Q', 'MaxItterQ', 'Reaxfflg', 'QEQInfo'), 
                              ('ReadVelocity', 'InitTemp', 'SeedTemp'),  
                              ("Thermo", "ThermoFreq"), 
                              ("Baro", "BaroFreq"), 
                              )
        pass
##################################################
    def set_path_wd(self, path:str):
        self.path_wd = Path.Path(path)
        return 0
##################################################
    def generate_unique_combis(self, length:int):
        """
        3体間や4体間の元素の組合せを作る。\n
        パラメータのAngleの項目やdump.bondからの結合の調査に利用。
        """
        # 重複を許可した組合せを生成
        type_set = set(self.atom["type"])
        combis = itertools.product(type_set, repeat=length)
        # 並びを反対にしたときに重複するならば削除
        unique_combis = []
        for combi in combis:
            reversed_combi = tuple(reversed(combi))
            if reversed_combi not in  unique_combis:
                unique_combis.append(combi)
        return unique_combis
##################################################
    def correct_dates(self):
        """_summary_
        cellsizeの最小値が負の場合、laichはerrorを発生させる。\n
        それ以外も解析をしやすくする目的でcellsizeの最小値を0に合わせる。\n
        それに応じてcellsizeの最大値や原子のxyz座標も調整する。\n
        inputfileを読み込み、cellsizeやatomsをinstanceが持っていれば実行可能。\n
        速度や力、マスク情報を持っていなければ全て0にして追加する。\n
        セルサイズ以上の値を持つ原子をセル範囲内に戻す。
        """
        minimums:np.array = self.cell[self.colm_cell[0]].copy()
        # positions
        for i, colum in enumerate(self.colm_atom[3:6]):
            self.atom[colum] -= minimums[i]
        # cellsize
        for colum in self.colm_cell:
            self.cell[colum] -= minimums
        # positions2
        minimums:np.array = self.cell[self.colm_cell[0]]
        maximums:np.array = self.cell[self.colm_cell[1]]
        positions = np.array([self.atom[colum] for colum in self.colm_atom[3:6]])
        for j, (minimum, maximum, qs) in enumerate(zip(minimums, maximums, positions)):
            qs_temp = np.array([0.]*qs.shape[0])
            for i, q in enumerate(qs):
                if minimum <= q < maximum:
                    qs_temp[i] = q
                elif q < minimum:
                    qs_temp[i] = q+maximum
                else:
                    qs_temp[i] = q-maximum
            self.atom[self.colm_atom[j+3]] = qs_temp.copy()
        # atom
        num_atom = self.atom[self.colm_atom[3]].shape[0]
        for colum in self.colm_atom:
            if not colum in self.atom.keys():
                self.atom[colum] = np.array([0.]*num_atom)
        del minimums, num_atom
        return 0
##################################################
    def correct_dtype(self):
        """_summary_
        atomsのid, type, maskはnp.int32型である必要がある。\n
        仕様上、読み込まれたままのデータは上記3つがnp.float64型なので、\n
        intに直す。
        """
        for colum in self.colm_atom[0:3]:
            self.atom[colum] = self.atom[colum].astype(np.int64)
        return 0
##################################################
    def reset_index(self):
        """
        atomのidデータをリセットする。
        """
        num_atom:int = self.atom["id"].shape[0]
        self.atom["id"] = np.array([i+1 for i in range(num_atom)])
        return 0
##################################################
    def make_mass(self):
        """_summary_
        self.massを作成する。
        self.atomsとself.elm_typeが存在する事が条件。
        よって、inputやdumoposと、para.rdを読み込んでいる必要がある。
        self.elm_typeを手動で{1,"C", 2:"H"}のように設定も可能。
        """
        type_set = set(self.atom[self.colm_atom[1]])
        self.mass = {
        self.colm_mass[0]: np.array([type_ for type_ in self.type_elm.keys() if type_ in type_set]), 
        self.colm_mass[1]: np.array([mass_dict[element] for type_, element in self.type_elm.items() if type_ in type_set])
        }
        return 0
##################################################
    def sort_by_id(self):
        """_summary_
        idの順番がバラバラになっている時、atomsをid順に直す。
        """
        sorted_index = np.argsort(self.atom['id'])
        temp_atoms = dict()
        for key, values in self.atom.items():
            temp_atoms[key] = [values[i] for i in sorted_index]
        self.atom = temp_atoms
        del sorted_index
        del temp_atoms
        return 0
##################################################
    def judge_numerical_conversion(self, string:str):
        """
        ある変数が数値変換可能かどうかを判定する。
        """
        try:
        # floatへの変換を試みる
            float(string)
            return True
        except ValueError:
            return False
##################################################
    def get_num_atom(self):
        """
       原子数を数えて返す 
        """
        return self.atom["id"].shape[0]
##################################################
    def make_periodic_atom(self) -> dict:
        """
        27方向の周期条件を考える為の原子情報の入ったdictを作成する
        """
        # 周期条件を考慮する27のセル方向を設定
        combis = np.array(list(itertools.product([1,0,-1], repeat=3)))
        
        # cellの長さを作成
        cell_max = self.cell["max"]
        cell_min = self.cell["min"]
        cell_len = cell_max - cell_min
        
        # 周期条件を考慮したatomを作成
        self.periodic_atom:dict = {key:np.array([]) for key in self.colm_atom}
        atom_origin = self.atom.copy()
        for combi in combis:
            add_len = cell_len * combi
            atom_temp_ = atom_origin.copy()
            for i, colm_q in enumerate(self.colm_atom[3:6]):
                atom_temp_[colm_q] = atom_origin[colm_q] + add_len[i]
            for colm in self.colm_atom:
                origin = self.periodic_atom[colm].copy()
                add = atom_temp_[colm].copy()
                self.periodic_atom[colm] = np.concatenate([origin, add], axis=0)
        for colum in self.colm_atom[0:3]:
            self.periodic_atom[colum] = self.periodic_atom[colum].astype(np.int64)
        return 0
##################################################
    def culc_moving_average(self, 
                            values:list, 
                            interbal:int=3, 
                            type_:str="SMA", 
                            decrement:float=1, 
                            ignore_part:str="nan") -> list:
        """
        type_:   default SMA
        SMA...単純移動平均
        WMA...加重移動平均
        いずれの場合も中央移動平均のみ求める
        values: 移動平均を求めたい数値群
        interbal: 移動平均を求める区間
        decrement: 減少率。大きい程、遠くの重みが小さくなる。
        ignore_part: nanならnp.nan、asなら元のリストの値
        """
        def get_weight(dist:int, type_:str=type_, decrement:float=decrement) -> float:
            if type_=="SMA":
                weight = 1.0
            elif type_=="WMA":
                weight = 1.0-(0.1*dist*decrement)
            return weight
        values:np.array = list(values)
        leng_values = len(values)
        
        if ignore_part=="nan":
            new_values:list = [np.nan for _ in range(leng_values)]
        elif ignore_part=="as":
            new_values:list = values.copy()
        
        if interbal%2==1:
            nums_ignore:int = interbal//2
            indices:list = range(leng_values)[nums_ignore:-nums_ignore]
            internal_indices:list = range(-nums_ignore, nums_ignore+1)
            for index in indices:
                value_temp:float = 0.0
                sum_weight:float = 0.0
                for internal_index in internal_indices:
                    weight:float = get_weight(dist=internal_index,
                                              type_=type_)
                    sum_weight += weight
                    value_temp += weight*values[index+internal_index]
                value_temp /= sum_weight
                new_values[index] = value_temp
            
        elif interbal%2==0:
            nums_ignore:int = interbal//2
            indices:list = range(leng_values)[nums_ignore:-nums_ignore]
            internal_indices:list = range(-nums_ignore+1, nums_ignore)
            for index in indices:
                value_temp:float = 0.0
                sum_weight:float = 0.0
                weight:float = get_weight(dist=nums_ignore,
                                          type_=type_)
                value_temp += 0.5*weight*values[index-nums_ignore]
                value_temp += 0.5*weight*values[index+nums_ignore]
                sum_weight += weight
                for internal_index in internal_indices:
                    weight:float = get_weight(dist=internal_index,
                                              type_=type_)
                    sum_weight += weight
                    value_temp += weight*values[index+internal_index]
                value_temp /= sum_weight
                new_values[index] = value_temp
        
        return new_values
##################################################
##################################################
##################################################