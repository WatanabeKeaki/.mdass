import numpy as np
import pathlib as Path
import pandas as pd
import re

from .cppfiles.judge_overlap import judge_overlap

##################################################
##################################################
class SelectAtom():
    def __init__(self) -> None:
        self.colm_cell:list = ["min", "max"]
        self.colm_mass:list = ["type", "mass"]
        self.colm_atom:list = ["id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q"]
        
        self.keys_int:list = ["id", "type", "mask"]
        pass
##################################################
    def select_dist(self, 
                    ids_remain:list, 
                    ids_judge:list, 
                    cutoff:float=2.0):
        """
        距離から原子を選択する
        """
        idxs_remain = [id_-1 for id_ in ids_remain]
        idxs_judge = [id_-1 for id_ in ids_judge]
        
        atom_remain = []
        atom_judge  = []
        for colm in ["id", "x", "y", "z"]:
            atom_remain.append(self.atom[colm][idxs_remain])
            atom_judge.append(self.atom[colm][idxs_judge])
        atom_remain = np.array(atom_remain).T.tolist()
        atom_judge = np.array(atom_judge).T.tolist()
        
        selected_ids = list(set(judge_overlap(atom_remain, atom_judge, cutoff)))
        self.selected_ids = list(map(int, selected_ids))
        return 0
##################################################
    def select_plane(self, 
                     norm_vect:list=[0.0,0.0,1.0], 
                     plane_point:list=[0.0,0.0,0.0], 
                     above_below:str="above",
                     standard:str="x"
                     ):
        assert len(norm_vect) == 3
        assert len(plane_point) == 3
                
        self.selected_ids:list = []
        
        dir_ind = {key: value for key, value in zip(["x","y","z"],[0,1,2])}
        standard_ind = dir_ind[standard]
        
        # assert int(norm_vect[standard_ind]) == 0.0
        
        # 切片dを決定
        d = 0
        norm_pp = list(zip(norm_vect, plane_point))
        for nrm, q in norm_pp:
            d += nrm*q
            
        # 計算に必要な法線ベクトル成分を残す
        rsdal_vect = norm_vect.copy()
        del rsdal_vect[standard_ind]
        
        # selected_idsにidを入れるか判別する
        ids:np.array = self.atom["id"]
        xyzs:list = []
        for colm_q in ["x","y","z"]:
            xyzs.append(self.atom[colm_q])
        xyzs:np.array = np.array(xyzs).T
        # xyzsは[[x0,y0,z0], [x1,y1,z1], ..., [xn,yn,zn]]の形式をとっている
        
        # 判定し、selected_idsをつくる
        for id_, xyz in zip(ids, xyzs):
            refer_xyz:np.array = xyz[standard_ind]
            check_xyz:np.array = xyz.copy()
            check_xyz = np.delete(check_xyz, standard_ind)
            # del check_xyz[standard_ind]
            refer_planepoint:float = (d-( xyz[0]*rsdal_vect[0] + xyz[1]*rsdal_vect[1] ))/norm_vect[standard_ind]
            
            if (above_below=="above" and refer_xyz>=refer_planepoint) or (above_below=="below" and refer_xyz<=refer_planepoint):
                self.selected_ids.append(id_)
        return 0
##################################################
    def select_type(self, select_types:list):
        """
        atom_typeにより選ぶ。
        引数のselect_typesはlistで与える
        例: [3]や[3, 5]
        self.selected_ids:listを返す。
        """
        selected_idx:list = []
        for select_type in select_types:
            selected_idx += np.where(self.atom["type"]==select_type)[0].tolist()
        self.selected_ids = [idx+1 for idx in selected_idx]
        return 0
##################################################
    def select_mask(self, select_masks:list):
        """
        maskにより選ぶ。
        引数のselect_masksはlistで与える
        例: [1]や[3, 2, ]
        self.selected_idsを返す。
        """
        selected_idx:list = []
        for select_mask in select_masks:
            selected_idx += np.where(self.atom["mask"]==select_mask)[0].tolist()
        self.selected_ids = [idx+1 for idx in selected_idx]
        return 0
##################################################
    def select_sin(self, 
                   amp_direction:str = "z", above_below:str = "above",
                   wavenum_i:int = 1, amplitude_i:float = None,
                   wavenum_j:int = 1, amplitude_j:float = None,
                   wave_intercept:float = None, 
                   print_info:bool = False
                   ):
        xyz = ["x","y","z"]
        dir_ind = {key: value for key, value in zip(xyz,[0,1,2])}
        standard_ind = dir_ind[amp_direction]
        del xyz[standard_ind]
        i, j = xyz
        i_direction = dir_ind[i]
        j_direction = dir_ind[j]
        
        atom_i:np.array = self.atom[i]
        atom_j:np.array = self.atom[j]
        atom_k:np.array = self.atom[amp_direction]
        
        len_i:np.array = self.cell["max"][i_direction]
        len_j:np.array = self.cell["max"][j_direction]
        len_k:np.array = self.cell["max"][standard_ind]
        
        if amplitude_i==None:
            amplitude_i = len_k/20
        if amplitude_j==None:
            amplitude_j = len_k/20
        if wave_intercept==None:
            wave_intercept = len_k/3
        
        if print_info:
            print(f"amplitude_max = {amplitude_i+amplitude_j}")
            print(f"wave_intercept = {wave_intercept}")
            print(f"minimum atlitude = {wave_intercept-(amplitude_i+amplitude_j)}")
            print(f"max atlitude = {wave_intercept+(amplitude_i+amplitude_j)}")
        
        def check_function(id_) -> float:
            value_i = amplitude_i*np.sin(2*np.pi*wavenum_i*atom_i[id_]/len_i)
            value_j = amplitude_j*np.sin(2*np.pi*wavenum_j*atom_j[id_]/len_j)
            check_k = value_i + value_j + wave_intercept
            return check_k
        
        self.selected_ids:list = []
        for idx, k in enumerate(atom_k):
            id_ = idx+1
            condition1 = (above_below=="above" and k>check_function(id_))
            condition2 = (above_below=="below" and k<check_function(id_))
            if condition1 or condition2:
                self.selected_ids.append(id_)
        
        return 0
##################################################
    def select_cos(self, 
                   amp_direction:str = "z", above_below:str = "above",
                   wavenum_i:int = 1, amplitude_i:float = None,
                   wavenum_j:int = 1, amplitude_j:float = None,
                   wave_intercept:float = None, 
                   print_info:bool = False
                   ):
        xyz = ["x","y","z"]
        dir_ind = {key: value for key, value in zip(xyz,[0,1,2])}
        standard_ind = dir_ind[amp_direction]
        del xyz[standard_ind]
        i, j = xyz
        i_direction = dir_ind[i]
        j_direction = dir_ind[j]
        
        atom_i:np.array = self.atom[i]
        atom_j:np.array = self.atom[j]
        atom_k:np.array = self.atom[amp_direction]
        
        len_i:np.array = self.cell["max"][i_direction]
        len_j:np.array = self.cell["max"][j_direction]
        len_k:np.array = self.cell["max"][standard_ind]
        
        if amplitude_i==None:
            amplitude_i = len_k/20
        if amplitude_j==None:
            amplitude_j = len_k/20
        if wave_intercept==None:
            wave_intercept = len_k/3
        
        if print_info:
            print(f"amplitude_max = {amplitude_i+amplitude_j}")
            print(f"wave_intercept = {wave_intercept}")
            print(f"minimum atlitude = {wave_intercept-(amplitude_i+amplitude_j)}")
            print(f"max atlitude = {wave_intercept+(amplitude_i+amplitude_j)}")
        
        def check_function(id_) -> float:
            value_i = amplitude_i*np.cos(2*np.pi*wavenum_i*atom_i[id_]/len_i)
            value_j = amplitude_j*np.cos(2*np.pi*wavenum_j*atom_j[id_]/len_j)
            check_k = value_i + value_j + wave_intercept
            return check_k
        
        self.selected_ids:list = []
        for idx, k in enumerate(atom_k):
            id_ = idx+1
            condition1 = (above_below=="above" and k>check_function(id_))
            condition2 = (above_below=="below" and k<check_function(id_))
            if condition1 or condition2:
                self.selected_ids.append(id_)
        
        return 0
##################################################
    def select_cylinder(self, radius:float=10, center:list=[0,0,0], c_axis:str="x"):
        """
        モデルを円柱状に選択する。円柱内部の原子idが保存される。
        c_axis_dirはc軸方向。円柱の高さ方向を指定する。"x"か"y", "z"を選ぶ。
        斜めなど角度のついた円柱を作りたい場合はrotation関数と組み合わせる。
        radiusは半径。初期値がでかすぎる、という場合は各自設定する。各軸cellsizeを超える直径は選択できない。
        centerは円の中心座標。周期条件は考慮されない。[x,y,z]のlistで指定する。
        reverseがFalseなら内部の原子が選択される。Trueなら内部が選択される。
        self.selected_idsを返す。
        """
        def culc_range(point0:np.array, point1:np.array) -> float:
            """
            平面座標の距離を測る。point0, point1は長さ2のlist
            """
            return np.linalg.norm(point0-point1)
        
        xyz:list = [None for _ in range(3)]
        for i, colm_q in enumerate(self.colm_atom[3:6]):
            xyz[i] = self.atom[colm_q].copy()
        xyz:np.array = np.array(xyz)
        
        # 判定に必要な要素のみを残す。c軸方向の情報を消す
        colm_xyz:list = list(self.colm_atom[3:6])
        direction_index = {colm_q: ind for ind, colm_q in enumerate(colm_xyz)}
        index = direction_index[c_axis]
        del center[index]
        
        colm_xyz.remove(c_axis)
        
        # 四角形で切り取る
        remaining_ids:list = []
        # condition0 = center[0]-radius<=xyz[]