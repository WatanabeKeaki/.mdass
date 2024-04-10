import sys
sys.path.append("./cppfiles/")
import dateset

import pickle
from multiprocessing import Pool
import numpy as np
from ase.build import molecule
import pathlib as Path
from time import time
import pandas as pd

np.load
class ImportAtom:
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
        
        self.cell
        self.masses
        self.atoms
        pass
##################################################
    def import_rd(self, name_input:str="input.rd") -> dict:
        """_summary_\n
        .rd形式のlaichで使われるinput fileを読み込む。
        Args:
            name_input (str, optional): _description_. Defaults to "input.rd".
        """
        # パスの作成
        path_input = self.path_wd.joinpath(name_input)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.split() for line in f.readlines(1000)]
        # 行ごとに解析
        for row, line in enumerate(lines):
            if not line:
                pass
            elif "#cellx" in line:
                cellx:np.array = np.array(line[1:]).astype(np.float64)
            elif "#celly" in line:
                celly:np.array = np.array(line[1:]).astype(np.float64)
            elif "#cellz" in line:
                cellz:np.array = np.array(line[1:]).astype(np.float64)
            elif "#masses" in line:
                masses_num:int = int(line[1])
                startrow_mass:int = row+1
            elif "#atoms" in line:
                atoms_num:int = int(line[1])
                startrow_atom:int = row+1
        # cellsize
        self.cell:dict = dict()
        dates_cell:np.array = np.array([cellx, celly, cellz]).T # .Tは転置
        for key, date in zip(self.colm_cell, dates_cell):
            self.cell[key] = np.array(date).astype(np.float64)
        # masses
        self.masses:dict = dict()
        dates_mass:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_mass, max_rows=masses_num, unpack=True)
        for key, date in zip(self.colm_mass, dates_mass):
            self.masses[key] = date
        # atoms
        self.atoms:dict = dict()
        dates_atom:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_atom, max_rows=atoms_num, unpack=True)
        for key, date in zip(self.colm_atom, dates_atom):
            self.atoms[key] = date
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
    def import_dump(self, name_input:str="dump.pos.0") -> dict:
        path_input = self.path_wd.joinpath(name_input)
        
        with open(path_input, "r") as f:
            lines = [line.split() for line in f.readlines(1000)]
        
        for row, line in enumerate(lines):
            if not line:
                pass
            elif "NUMBER OF ATOMS" in " ".join(line):
                atoms_num:int = int(lines[row+1][0])
            elif "BOX BOUNDS" in " ".join(line):
                cells_num = 3
                startrow_cell = row+1
            elif "ATOMS id type" in " ".join(line):
                startrow_atom = row+1
                colm_atom = line[2:]
                break
        # cellsize
        self.cell:dict = dict()
        dates_cell:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_cell, max_rows=cells_num, unpack=True)
        for key, date in zip(self.colm_cell, dates_cell):
            self.cell[key] = date
        # atoms
        self.atoms:dict = dict()
        dates_atom:np.array = np.loadtxt(fname=str(path_input), skiprows=startrow_atom, max_rows=atoms_num, unpack=True)
        sorted_indices = np.argsort(dates_atom[0])
        dates_atom = dates_atom[:, sorted_indices]
        for key, date in zip(colm_atom, dates_atom):
            self.atoms[key] = date
        if self.colm_atom[2] not in colm_atom:
            masks = np.array([0*atoms_num])
            self.atoms[self.colm_atom[2]] = masks
            
        # correction
        self.correct_dates()
        self.correct_dtype()
        # masses
        self.make_masses()
        return 0
##################################################
##################################################
class ImportPara:
    def __init__(self) -> None:
        pass
    
    def import_para(self, name_input:str="para.rd") -> dict:
        """
        self.para_dictとself.parahead_dictを出力する
        
            para_dict={
                atom:{(C:str,):[100, 0.1, ...], (H,):[120, 0.2, ...], 
                bond:{(1:int, 2:int):[100, 0.1, ...]}
                ...}
            の形式でパラメータを格納する
            parahead_dict:
                {
                atom:"Reactive MD-force field for C/H/O/N/Si published version developed by yang wang\n...",
                bond:" 15       ! Nr of bonds; Edis1;LPpen;n.u.;pbe1;pbo5;13corr;pbo6\n"pbe2;...
                ...
                }
            の形式でheaderを格納する
        
        """
        path_input = self.path_wd.joinpath(name_input)
        # para.rdを読み込む為に使用
        checkword_list = ['Nr of atoms', 'Nr of bonds','Nr of off-diagonal', 
                          'Nr of angles', 'Nr of torsions', 'Nr of hydrogen bonds']
        #パラメータを格納するdictに利用
        parakey_list:list = ["atom", "bond", "off-diagonal", 
                             "angle", "torsion", "hydrogen-bond"]
        #fileの読み込み
        with open(path_input, "r") as f:
            lines = [line.split() for line in f.readlines()]
        
        # AtomやBondの項目が何行目から始まるか確認し、inditem_listに格納する
        inditem_dict = dict(zip(parakey_list, [None]*len(parakey_list)))
        # AtomやBondの中にいくつの組み合わせ等があるか格納する
        numitem_dict = dict(zip(parakey_list, [None]*len(parakey_list)))
        
        # inditem_dictとnumitem_dictを決定
        for para_index, line in enumerate(lines):
            for checkword, key_str in zip(checkword_list, parakey_list):
                if checkword in " ".join(line):
                    inditem_dict[key_str] = (para_index)
                    numitem_dict[key_str] = int(line[0])
        
        # para_dictに保存されないatomやbond, offdiagonalのheader情報を保存する
        self.parahead_dict = dict(zip(parakey_list, [None for _ in parakey_list]))
        Num_Head_tuple = (4, 2, 1, 1, 1, 1)
        startind_list = list(inditem_dict.values())
        for key_index, key_str in enumerate(parakey_list):
            if key_index == 0:
                head_list = lines[0:startind_list[key_index]+Num_Head_tuple[key_index]].copy()
            else:
                head_list = lines[startind_list[key_index]:startind_list[key_index]+Num_Head_tuple[key_index]].copy()
            para_heads = "" # 初期化
            for head in head_list:
                para_heads += " ".join(head)+"\n"
            # print(para_heads) # テスト用
            self.parahead_dict[key_str] = para_heads
        # print(self.parahead_dict)
        
        # パラメータで変更する可能性がある情報を保存する。具体的なparamaterを保存
        self.para_dict = dict(zip(parakey_list, [{} for _ in parakey_list]))
        num_colum_tuple = (1, 2, 2, 3, 4, 3)
        def store_para(key_index, key_str):
            ind_start = inditem_dict[parakey_list[key_index]]
            num_spcs = numitem_dict[parakey_list[key_index]]
            for spc in range(1, num_spcs+1):
                for ind in range(Num_Head_tuple[key_index]):
                    num_skip:int = num_colum_tuple[key_index]
                    line = lines[ind_start+spc*Num_Head_tuple[key_index]+ind]
                    if ind == 0:
                        elm_tuple:tuple = tuple(line[0:num_skip])
                        self.para_dict[key_str][elm_tuple] = line[num_skip:]
                    else:
                        self.para_dict[key_str][elm_tuple] += line
            return 0
        for key_index, key_str in enumerate(parakey_list):
            store_para(key_index, key_str)
        
        #atomsのtypeとparaの情報を一致させる
        elms:list = list(self.para_dict[parakey_list[0]].keys())
        def get_elmtype() -> dict:
            for i, elm_ in enumerate(elms):
                elms[i] = elm_[0]
            type_list = [i+1 for i in range(len(elms))]
            self.elm_type = dict(zip(elms, type_list))
            self.type_elm = dict(zip(type_list, elms))
            return self.elm_type, self.type_elm
        self.elm_type:dict
        self.type_elm:dict
        get_elmtype()
        return 0
##################################################
##################################################
class ImportBond:
    def __init__(self) -> None:
        pass
##################################################
    def import_bond(self, name_input:str="dump.bond.0"):
        times = []
        times.append(time())
        path_input = self.path_wd.joinpath(name_input)
        
        with open(path_input, "r") as f:
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
        
        #律速
        times.append(time())
        df_bond:pd.DataFrame = pd.read_csv(path_input, skiprows=skiprows, 
                                           sep=" ", header=None, 
                                           names=[0, 1, 2, 3, 4], engine="c")
        df_bond:pd.DataFrame = df_bond.replace('Atom', -1).astype("float32")
        array_bond:np.array = df_bond.values
        del df_bond
        # array_bond = np.array(get_bond(str(path_input)))
        
        # bondorderを抽出
        # times.append(time())
        bondorders:np.array = array_bond[:, 1:5]
        bondids = array_bond[:, 0]
        
        # 欠損値を含む行を選択するマスクを作成する,bool型をvalueに持つnp.arrayが作成される
        mask_nan = np.any(np.isnan(array_bond), axis=1)
        rownums_mainhead = np.where(mask_nan)[0]
        
        # 欠損値を含む行の2番目と3番目の列を選択し、整数型に変換する
        mainid_bondnum = array_bond[mask_nan, 1:3].astype(int).T
        mainids = mainid_bondnum[0]
        bondnums = mainid_bondnum[1]
        
        # 第二律速
        # 以下2つのarrayは一対一対応でなければならない
        # times.append(time())
        self.bondorder_list = [None]*num_atoms
        self.bondorder_connect_list = [None]*num_atoms
        for rownum, mainid, bondnum in zip(rownums_mainhead, mainids, bondnums):
            self.bondorder_list[mainid-1] = bondorders[rownum+1: rownum+bondnum+1]
            self.bondorder_connect_list[mainid-1] = bondids[rownum+1: rownum+bondnum+1]
        
        # times.append(time())
        # for i in range(len(times)-1):
        #     print(times[i+1]-times[i])
        
        
        return 0