import numpy as np
import pathlib as Path
import pandas as pd
import re
from time import time

##################################################
##################################################
class OperateAtom():
    def __init__(self) -> None:
        self.atom: dict
        
        self.colm_atom:list = ("id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q")
        self.colm_cell:list = ("min", "max")
        self.colm_mass:list = ("type", "mass")
        pass
##################################################
    def remain_atom(self):
        """
        必要なデータ:
        self.atom... import_input()等から得られる。
        selected_ids:list... select_sin(), select_plane()等から得られる。
        説明：
        remove_atoms()では、selected_ids内のidに属する原子のみをcell内から除外。
        使用例: 摺動モデルを作る際、選択した原子のみ除外
        """
        selected_indices = np.array(self.selected_ids)-1
        atom_temp_:dict = dict()
        for key, value in self.atom.copy().items():
            atom_temp_[key] = value[selected_indices]
        self.atom:dict = atom_temp_.copy()
        self.reset_index()
        del atom_temp_
        return 0
##################################################
    def remove_atom(self):
        """
        必要なデータ:
        self.atom... import_input()等から得られる。
        selected_ids:list... select_sin(), select_plane()等から得られる。
        説明：
        remove_atoms()では、selected_ids内のidに属する原子のみをcell内から除外。
        使用例: 摺動モデルを作る際、選択した原子のみ除外
        """
        selected_indices = np.array(self.selected_ids)-1
        all_indices = self.atom["id"]-1
        remaining_indices = np.setdiff1d(all_indices, selected_indices)
        atom_temp_:dict = dict()
        for key, value in self.atom.copy().items():
            atom_temp_[key] = value[remaining_indices]
        self.atom:dict = atom_temp_.copy()
        self.reset_index()
        del atom_temp_
        return 0
##################################################
    def change_mask(self, mask:int):
        """
        必要なデータ:
        self.atoms_df... import_input()等から得られる。
        selected_ids:list... select_sin(), select_plane()等から得られる。
        説明：
        change_mask()では、selected_ids内のidに属する原子のmaskを変更可能。
        使用例: 注意して観測したい原子のみmaskを指定する
        """
        selected_idxs = [id_-1 for id_ in self.selected_ids]
        for idx in selected_idxs:
            self.atom["mask"][idx] = mask
        return 0
##################################################
    def change_type(self, type_:int):
        """
        必要なデータ:
        self.atoms_df... import_input()等から得られる。
        selected_ids:list... select_sin(), select_plane()等から得られる。
        説明：
        change_type()では、selected_ids内のidに属する原子のtypeを変更可能。
        使用例: FeバルクのFe原子の一部をTi等と置換し、置換型固溶体を作成する
        """
        selected_idxs = [id_-1 for id_ in self.selected_ids]
        for idx in selected_idxs:
            self.atom["type"][idx] = type_
        return 0
##################################################
    def expand_cell(self, replicate:list = [1, 1, 1]):
        """
        replicate               : セルをreplicateにしたがって拡大する。[1,2,2]など３要素のリストで与える
        """
        len_x, len_y, len_z = self.cell["max"]
        reps_x, reps_y, reps_z = replicate
        
        # cellsize
        self.cell["max"] *= replicate
        
        # atom
        atom_new:dict = dict()
        for colm in self.colm_atom:
            atom_new[colm] = np.array([])
        
        # expand
        atom_origin:dict = self.atom.copy()
        for rep_x in range(reps_x):
            for rep_y in range(reps_y):
                for rep_z in range(reps_z):
                    atom_temp_ = atom_origin.copy()
                    
                    addvalue_x = rep_x*len_x
                    atom_temp_["x"] = atom_origin["x"] + addvalue_x
                    
                    addvalue_y = rep_y*len_y
                    atom_temp_["y"] = atom_origin["y"] + addvalue_y
                    
                    addvalue_z = rep_z*len_z
                    atom_temp_["z"] = atom_origin["z"] + addvalue_z
                    
                    for colm in self.colm_atom:
                        origin = atom_new[colm]
                        add = atom_temp_[colm]
                        atom_new[colm] = np.concatenate((origin, add))
        self.atom = atom_new.copy()
        
        # correction
        self.reset_index()
        self.correct_dtype()
        return 0
##################################################
    def extend_cell(self, rates:list=[1.0,1.0,1.0]):
        """
        原子数等を変えずに、セルを引き延ばす。
        密度やセルサイズの調整のために用いる
        rates: x,y,z方向にどれくらい拡張するかを記入。
        """
        self.cell["max"] *= rates
        
        for colm_q, rate in zip(self.colm_atom[3:6], rates):
            self.atom[colm_q] *= rate
        return 0
##################################################
    def culc_density(self) -> float:
        """
        計算した密度を返す
        """
        # セルの体積を計算
        volume = self.cell["max"].prod()*(10**-24)
        
        # セル内の原子の重さの和をとる
        sum_ = 0.0
        types = self.mass["type"]
        masses = self.mass["mass"]
        for type_, mass in zip(types, masses):
            indices = np.where(self.atom["type"]==type_)[0]
            num_atom = indices.shape[0]
            sum_ += num_atom*mass/(6.02214076*10**23)
        density = sum_ / volume
        print(density)
        return density
##################################################
    def marge_atom(self, instances:list):
        """
        2つ以上のatomsやmassesの情報を結合する。
        cellsizeはinstancesの一番目が採用される
        massesのatomtypeは統一されていなければならない。(Oなら3, Siなら4など)
        instancesはinstanceのリスト。
        """
        # atom
        for colm in self.colm_atom:
            values_:list = [None for _ in range(len(instances))]
            for i, instance in enumerate(instances):
                values_[i] = instance.atom[colm]
            self.atom[colm] = np.concatenate(values_)
        # mass
        self.make_mass()
        # cell
        self.cell = instances[0].cell
        # correction
        self.reset_index()
        self.correct_dtype()
        self.correct_dates()
        return 0
##################################################
    def reverse_atom(self, direction:str="z"):
        """
        direction方向について反転する。
        cellsizeは保たれる
        direction: "x","y","z"
        のいずれか。
        """
        cellmax_dict = {key: value for key, value in zip(["x","y","z"],self.cell["max"])}
        self.atom[direction] = -1*self.atom[direction].copy()
        self.atom[direction] = self.atom[direction].copy() + cellmax_dict[direction]
        return 0
##################################################
    def add_space(self, amount:float, direction:str="z"):
        """
        axisで指定した方向に、amountで指定した長さのスペースを追加する
        """
        
        dir_ind = {key: value for key, value in zip(["x","y","z"],[0,1,2])}
        ind = dir_ind[direction]
        self.cell["max"][ind] += amount
        return 0
##################################################
    def rotation(self, center:list=[0,0,0], eulerian_angles:list=[0,0,0], margin:float = 2):
        """
        center          : 回転中心の座標, 3次元のみ。xyz座標
        eulerian_angles : オイラー角。ラジアン角を使用。
        margin          : 回転後に原子が周期境界を超える事を防ぐ為に設定
        """
        center = np.array(center)
        
        fai:float = eulerian_angles[0]
        theta:float = eulerian_angles[1]
        psi:float = eulerian_angles[2]
        
        rot_mat_x = np.array([[1,           0,            0],
                              [0, np.cos(fai), -np.sin(fai)],
                              [0, np.sin(fai),  np.cos(fai)]])
        rot_mat_y = np.array([[ np.cos(theta), 0, np.sin(theta)],
                              [             0, 1,             0],
                              [-np.sin(theta), 0, np.cos(theta)]])
        rot_mat_z = np.array([[np.cos(psi), -np.sin(psi), 0],
                              [np.sin(psi),  np.cos(psi), 0],
                              [          0,            0, 1]])
        
        rot_mat = rot_mat_x @ rot_mat_y @ rot_mat_z
        
        xyzs = []
        for colm_q in self.colm_atom[3:6]:
            qs = self.atom[colm_q]
            xyzs.append(qs)
        xyzs = np.array(xyzs).T
        
        num_atom = self.atom["id"].shape[0]
        
        xyzs_rot:list = [0]*num_atom
        for idx, xyz in enumerate(xyzs):
            xyzs_rot[idx] = rot_mat @ xyz
        xyzs_rot = np.array(xyzs_rot).T
        
        xyz_min = []
        for qs in xyzs_rot:
            xyz_min.append(np.min(qs))
        xyzs_rot = np.array(xyzs_rot).T
        xyz_min = np.array(xyz_min)
        
        xyzs_new = [0]*num_atom
        for idx, xyz_rot in enumerate(xyzs_rot):
            xyz_new = xyz_rot-xyz_min+margin
            xyzs_new[idx] = xyz_new.copy()
        xyzs_new = np.array(xyzs_new).T
        xyz_max = []
        for qs, colm_q in zip(xyzs_new, self.colm_atom[3:6]):
            xyz_max.append(np.max(qs))
            self.atom[colm_q] = qs
        
        self.cell["min"] = np.array([0., 0., 0.])
        self.cell["max"] = np.array(xyz_max)+margin
        
        return 0
##################################################
    def embed_atom(self, 
                   main, 
                   embed, 
                   start_setting:str="origin", 
                   coordinate:list=[0,0,0], 
                   cutoff:float=2.0):
        """
        バルク中に粒子を埋め込む際に使用する
        cellsizeは変化しない
        main_struct     バルクのHmFrame
        embedded_struct 粒子のHmFrame
        starting_point  粒子構造ファイルのどこを始点にするか
                        origin: 原点
                        center: 中央
        coordinate      始点の座標
        """
        instances:list = [main, embed]
        
        num_atom_main:int = main.get_num_atom()
        num_atom_embed:int = embed.get_num_atom()
        ids_remain:list = [i+1 for i in range(num_atom_main, num_atom_main+num_atom_embed)]
        
        #削除する原子の範囲を設定
        coordinate = np.array(coordinate)
        cell_emb = embed.cell["max"]
        if start_setting=="origin":
            start_point = coordinate
            finish_point = cell_emb
        elif start_setting == "center":
            start_point = coordinate-(cell_emb*0.5)
            finish_point = coordinate+(cell_emb*0.5)
        
        atom_emb = embed.atom.copy()
        for colm_q, sp in zip(self.colm_atom[3:6], start_point):
            atom_emb[colm_q] += sp
        
        main.make_periodic_atom()
        
        ids:np.array = main.periodic_atom["id"]
        xs_mp:np.array = main.periodic_atom["x"]
        ys_mp:np.array = main.periodic_atom["y"]
        zs_mp:np.array = main.periodic_atom["z"]
        xyzs_mp:np.array = np.array([xs_mp, ys_mp, zs_mp]).T
        
        min_value:np.array = start_point - cutoff
        max_value:np.array = finish_point + cutoff
        cons = np.all((xyzs_mp >= min_value) & (xyzs_mp <= max_value), axis=1)
        ids_judge:list = ids[cons].tolist()
        
        self.marge_atom(instances=instances)
        self.select_dist(ids_remain=ids_remain,
                         ids_judge=ids_judge,
                         cutoff=cutoff)
        self.remove_atom()
        
        self.selected_ids = []
        return 0
##################################################
    def alternate_axis(self, axis0:str, axis1:str):
        """
        軸を入れ替えるときに使う
        position, valocity, force, cellsizeを入れ替える
        """
        dir_ind = {key: value for key, value in zip(["x","y","z"],[0,1,2])}
        ind_0 = dir_ind[axis0]
        temp_cell_min = self.cell["min"][ind_0]
        temp_cell_max = self.cell["max"][ind_0]
        
        temp_qs = self.atom[axis0].copy()
        temp_vel = self.atom[f"v{axis0}"].copy()
        temp_force = self.atom[f"f{axis0}"].copy()
        
        self.atom[axis0] = self.atom[axis1].copy()
        self.atom[f"v{axis0}"] = self.atom[f"v{axis1}"].copy()
        self.atom[f"f{axis0}"] = self.atom[f"f{axis1}"].copy()
        ind_1 = dir_ind[axis1]
        self.cell["min"][ind_0] = self.cell["min"][ind_1]
        self.cell["max"][ind_0] = self.cell["max"][ind_1]
        
        self.atom[axis1] = temp_qs
        self.atom[f"v{axis1}"] = temp_vel
        self.atom[f"f{axis1}"] = temp_force
        self.cell["min"][ind_1] = temp_cell_min
        self.cell["max"][ind_1] = temp_cell_max
        return 0
##################################################
##################################################
##################################################
class Operate(OperateAtom):
    def __init__(self) -> None:
        super().__init__()
        pass