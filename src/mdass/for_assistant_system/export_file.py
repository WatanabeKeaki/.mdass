import pandas as pd
import numpy as np
import pathlib as Path
from collections import Counter
from datetime import datetime
import csv
from time import time

class ExportAtom:
    def __init__(self) -> None:
        self.colm_cell:list = ("min", "max")
        self.colm_mass:list = ("type", "mass")
        self.colm_atom:list = ("id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q")
        pass
##################################################
    def export_xyz(self, fname:str, type_to_element:bool=True, 
                   id_:bool=True, velocity_:bool=True, 
                   mask_:bool=True, force_:bool=False, 
                   charge_:bool=False, velo_mag_:bool=False, 
                   packmol:bool = False):
        """_summary_
        xyzの形式で出力する。\n
        elementがTrueなら、typeが数字でなく元素記号になる
        Args:
            fname (str): _description_
        """
        # 出力先パス
        fpath = self.path_wd.joinpath(fname)
        # packmol=Trueなら、自動的に書き込まれる要素は決定される
        if packmol:
            type_to_element = True
            id_,velocity_,mask_,force_,charge_,velo_mag_=False,False,False,False,False,False
        # 出力項目の決定
        columns:list = []
        properties_selected:list = []
        properties_all:tuple = ('id:I:1',     'species:S:1', 'pos:R:3',
                                 'velo:R:3',   'mask:R:1',    'force:R:3',
                                 'charge:R:1', 'velo_mag:R:1')
        if id_:
            columns += self.colm_atom[0:1] # ["id"]
            properties_selected += properties_all[0:1]
        columns += self.colm_atom[1:2]+self.colm_atom[3:6] # ["type", "x", "y", "z"]
        properties_selected += properties_all[1:3]
        if velocity_:
            columns += self.colm_atom[6:9] # ["vx", "vy", "vz"]
            properties_selected += properties_all[3:4]
        if mask_:
            columns += self.colm_atom[2:3] # ["mask"]
            properties_selected += properties_all[4:5]
        if force_:
            columns += self.colm_atom[9:12] # ["fx", "fy", "fz"]
            properties_selected += properties_all[5:6]
        if charge_:
            columns += self.colm_atom[12:13] # ["q"]
            properties_selected += properties_all[6:7]
        if velo_mag_:
            pass
        # make header part
        headers = []
        # num_atom
        num_atom:int = self.atom[self.colm_atom[1]].shape[0]
        headers.append(f"{num_atom}\n")
        # cell
        len_x = self.cell[self.colm_cell[1]][0]
        len_y = self.cell[self.colm_cell[1]][1]
        len_z = self.cell[self.colm_cell[1]][2]
        lattice:str = f'Lattice="{len_x:.4f} 0.0 0.0 0.0 {len_y:.4f} 0.0 0.0 0.0 {len_z:.4f}"'
        headers.append(f"{lattice}\t")
        # properties
        properties_selected:str = "Properties="+":".join(properties_selected)
        headers.append(f"{properties_selected}\n")
        # make csv part
        atom_temp:dict = {colum:self.atom[colum].copy() for colum in columns}
        if type_to_element:
            types_temp = atom_temp[self.colm_atom[1]]
            elements_temp = np.array([self.type_elm[type_] for type_ in types_temp])
            atom_temp[self.colm_atom[1]] = elements_temp.copy()
        out_df = pd.DataFrame(atom_temp)
        
        # write
        with open(file=fpath, mode="w") as f:
            for header in headers:
                f.write(header)
        out_df.to_csv(path_or_buf=fpath, sep="\t", header=None, 
                      index=None, mode="a")
        del fpath, columns, properties_selected, num_atom, len_x, len_y, len_z , lattice, headers, atom_temp, out_df
        return 0
##################################################
    def export_rd(self, fname:str="input.rd", options:list[str]=[], velocity_=True):
        """_summary_
        input.rdを作成する。

        Args:
            fname (str):\n
            出力するファイルの名前。\n
            options (list[str]):\n
            fixやmoveなどを書き加える。リスト内のstringsを直接inputに反映する。\n
            頭に#を書き忘れないよう注意。
        """
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        # 出力項目の決定
        columns = self.colm_atom[0:6]
        if velocity_:
            columns += self.colm_atom[6:9]
        # make header part
        headers = []
        # cell
        for i, q in enumerate(["x","y","z"]):
            cellmin = self.cell[self.colm_cell[0]][i]
            cellmax = self.cell[self.colm_cell[1]][i]
            headers.append(f"#cell{q}\t{cellmin:.4f}\t{cellmax:.4f}\n")
        # mass
        num_mass:int = len(self.mass[self.colm_mass[0]])
        headers.append(f"\n#masses {num_mass}\n")
        for i in range(num_mass):
            type_ = self.mass[self.colm_mass[0]][i]
            mass_ = self.mass[self.colm_mass[1]][i]
            headers.append(f"{type_}\t{mass_}\n")
        # option
        headers.append(f"\n")
        for option in options:
            headers.append(f"{option}\n")
        # atom
        num_atom = self.atom[self.colm_atom[0]].shape[0]
        headers.append(f"\n#atoms {num_atom}\n")
        # make csv part
        atom_temp:dict = {colum:self.atom[colum].copy() for colum in columns}
        out_df = pd.DataFrame(atom_temp)
        
        # write
        with open(fpath, "w") as f:
            for header in headers:
                f.write(header)
        out_df.to_csv(path_or_buf=fpath, sep="\t", header=None, 
                      index=None, mode="a")
        del fname, options, velocity_, fpath, columns, cellmin, cellmax, num_mass, type_, mass_, num_atom, atom_temp, out_df
        return 0
##################################################
    def export_POSCAR(self, fname:str="POSCAR", header:str="None"):
        """
        headerは1行目のコメントラインを指定。デフォルトはNone。
        """
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        lines:list = list()
        #header
        lines.append(header)
        #cell
        cell_max = self.cell[self.colm_cell[1]]
        cellx = cell_max[0]
        celly = cell_max[1]
        cellz = cell_max[2]
        lines.append("1.0")
        lines.append(f"{cellx} 0.0 0.0")
        lines.append(f"0.0 {celly} 0.0")
        lines.append(f"0.0 0.0 {cellz}")
        #element, num
        types_temp:np.array = self.atom[self.colm_atom[1]].copy()
        elements_temp:list = [self.type_elm[type_] for type_ in types_temp]
        elm_set = sorted(set(elements_temp))
        counts_elm = Counter(elements_temp)
        self.elm_num = {elm:str(counts_elm[elm]) for elm in elm_set}
        lines.append(" ".join(self.elm_num.keys()))
        lines.append(" ".join(self.elm_num.values()))
        #pos
        # #debug
        # self.atom["mask"] = self.atom["type"].copy()
        #end
        lines.append("Cartesian")
        atom_temp:np.array = np.array([self.atom[colum].copy() for colum in self.colm_atom[3:6]])
        standard_row = self.atom[self.colm_atom[1]].copy()
        sort_order = np.argsort(standard_row)
        atom_temp = atom_temp[:, sort_order]
        # write
        with open(fpath, "w") as f:
            for line in lines:
                f.write(f"{line}\n")
            np.savetxt(f, atom_temp.T, delimiter=' ', fmt='%.10f')
        return 0
##################################################
    def export_car(self, fname:str, PBC=True):
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        
        now = datetime.now()
        output_date = now.strftime("%a %b %d %H:%M:%S %Y\n")
        headers = [
            "!BIOSYM archive 3\n",
            "", # PBC
            "Materials Studio Generated CAR File\n",
            f"!DATE {output_date}"
        ]
        if PBC:
            headers[1] = "PBC=ON\n"
            cell_max = self.cell[self.colm_cell[1]]
            cellx = cell_max[0]
            celly = cell_max[1]
            cellz = cell_max[2]
            headers.append(
                f"PBC {cellx:8.4f} {celly:8.4f} {cellz:8.4f}    90.0000   90.0000   90.0000 (P1)\n")
        else:
            headers[1] = "PBC=OFF\n"
        
        
        types_temp:np.array = self.atom[self.colm_atom[1]].copy()
        elements_temp:np.array = np.array([self.type_elm[type_] for type_ in types_temp])
        counts_elm = dict(Counter(elements_temp))
        
        elememts_ind:list = []
        xyz:list = [
            [], # x
            [], # y
            []  # z
        ]
        unknown_strs:list = [] # carファイル特有のよくわからん文字列を格納
        elements:list = []
        unknown_nums:list = []
        
        for elem, num in counts_elm.items():
            elememts_ind += [f"{elem+str(ind+1):6}" for ind in range(num)]
            indices = np.where(elements_temp == elem)[0]
            xyz_temp:list = [list(self.atom[q][indices]) for q in self.colm_atom[3:6]]
            for i, qs_temp in enumerate(xyz_temp):
                xyz[i] += [f" {q:13.9f} " for q in qs_temp]
            unknown_strs += ["XXXX 1      xx      "]*num
            elements += [f"{elem:2}"]*num
            unknown_nums += ["  0.000"]*num
        xs = xyz[0]
        ys = xyz[1]
        zs = xyz[2]
        combined_lists:list = [
            elememts_ind, 
            xs, ys, zs, 
            unknown_strs, 
            elements, 
            unknown_nums
        ]
        transposed_lists = list(map(list, zip(*combined_lists)))
        lines = ["".join(line)+"\n" for line in transposed_lists]
        # write
        with open(file=fpath, mode="w") as f:
            for header in headers:
                f.write(header)
            for line in lines:
                f.write(line)
        return 0
##################################################
##################################################
class ExportPara:
    def __init__(self) -> None:
        pass
##################################################
    def export_para(self, fname:str="para.rd"):
        path_exportfile = self.path_wd.joinpath(fname)
        
        pnum_per_index_list = [
                    [9, 8, 8, 8], #atom
                    [10, 8], #bond
                    [8], #offdiagonal
                    [10], #angle
                    [11], #torsion
                    [7], #hydrogenbond
                              ]
        
        def decide_space(value:str, leng:int = 10) -> str:
            spacenum = leng - len(value)
            return " "*spacenum
        
        def decide_atomspace(elm):
            if len(elm) == 1 and len(elm[0]) == 1:
                return " "
            else:
                return ""
        
        para_lines:list = ["" for _ in pnum_per_index_list]
        for i, item_dict in enumerate(self.para_dict.values()):
            item_dict:dict
            
            for elm, elm_para_list in item_dict.items():
                elm:tuple
                elm_para_list:list
                
                atom_space = decide_atomspace(elm)
                para_lines[i] += "  "+"  ".join(elm)+atom_space
                init_space = " "*len("  "+"  ".join(elm)+atom_space)
                
                for j,para in enumerate(elm_para_list):
                    if j%8 == 0 and j != 0:
                        para_lines[i] += "\n"+init_space
                    
                    spaces = decide_space(para)
                    para_lines[i] += spaces+para
                    
                    if j+1 == len(elm_para_list):
                        para_lines[i] += "\n"
        
        write_lines:list = []
        for para_head, item_dict in zip(self.parahead_dict.values(), para_lines):
            write_lines.append(para_head)
            write_lines.append(item_dict)
        
        with open(path_exportfile, "w") as f:
            for line in write_lines:
                f.write(line)
                
        return 0
##################################################
    def export_POTCAR(self, 
                      fname="POTCAR", rootpath:str="/nfshome15/POTLIST/", 
                      functional:str="PBE", number = "54", 
                      options:dict = {"None": None}
                      ):
        """
        fnameは出力ファイル名;\n
        default: "POTCAR"\n
        \n
        rootpathはPOTCARを格納するディレクトリ;\n
        default: "/nfshome15/POTLIST/"\n
        \n
        functionalはどの汎関数を使うかを選択する;\n
        functional = [ "PBE" | "GGA" | "LDA" ]\n
        default: "PBE"\n
        \n
        numberはPBE.52, PBE.54, PBE.64のうちどの番号を使うか指定\n
        default: 54\n
        \n
        optionsはdict形式でどのようなpotentialを使うか指定する\n
        format: {element0: GW, element1: sv, element2: d, ...}\n
        example: {H: h, O: h}\n
        default: {"None": None}\n
        ※ hは短い結合を有す分子が存在する場合に使用が推奨される。\n
        どのようなpotentialを使うかは、対象の物質を扱った第一原理計算を用いる論文を調べ、\n
        計算の対象としている価電子数を見るのが最も確実。\n
        詳しくは、https://www.vasp.at/wiki/index.php/Available_PAW_potentials\n
        """
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        lines = []
        name_functional:str = f"potpaw_{functional}.{number}"
        path_POTCARdir:Path.Path = Path.Path(rootpath).joinpath(name_functional)
        for element in self.elm_num.keys():
            if element in options.keys():
                type_potential:str = options[element]
                potential:str = f"{element}_{type_potential}"
            else:
                potential:str = element
            path_potential:Path.Path = path_POTCARdir.joinpath(potential, "POTCAR")
            with open(path_potential, "r") as f:
                lines += f.readlines()
        with open(fpath, "w") as f:
            for line in lines:
                f.write(line)
        del fpath, lines, name_functional, path_POTCARdir
        return 0
##################################################
##################################################
class ExportMisc:
    def __init__(self) -> None:
        self.keywords_config_laich = (
            
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
        
        self.keywords_vasp = (
            ("SYSTEM")
        )
        pass
##################################################
    def export_config_laich(self, fname="config.rd"):
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        lines:list = []
        maxlen_key = max([len(key) for key in self.incar.keys()])
        for keywords in self.keywords_config_laich:
            for key in keywords:
                if key in self.config_laich.keys():
                    value = self.config_laich[key]
                    keylen = len(key)
                    spaces = " "*(maxlen_key-keylen+1)
                    lines.append(f"{key}{spaces}{value}\n")
            lines.append("\n")
        with open(fpath, "w") as f:
            for line in lines:
                f.write(line)
        del fpath, lines, value, keylen, spaces
        return 0
##################################################
    def export_KPOINTS(self, fname:str="KPOINTS", option:str="Gamma"):
        """
        今後オプションで色々指定できるようにする
        """
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        if option=="Gamma":
            lines:str="Gamma\n0\nGamma\n1\t1\t1"
        with open(fpath, "w") as f:
            for line in lines:
                f.write(line)
        return 0
##################################################
    def export_INCAR(self, fname:str="INCAR"):
        # 出力先パス
        fpath:Path.Path = self.path_wd.joinpath(fname)
        lines:list = []
        maxlen_key = max([len(key) for key in self.incar.keys()])
        for key, value in self.incar.items():
            keylen = len(key)
            spaces = " "*(maxlen_key-keylen) #タグ(key)の中で最長は "LPHON READ FORCE CONSTANTS"
            lines.append(f"{key}{spaces} = {value}\n")
        with open(fpath, "w") as f:
            for line in lines:
                f.write(line)
        del fpath, lines, value, keylen, spaces
        return 0
##################################################
    def export_jobscript_masamune(self, 
                                  task:str="MD", 
                                  nodes:int=1, 
                                  walltime:float=None, 
                                  que:str="P_016", 
                                  jobname:str="watanabe", 
                                  emailaddress:str="watanabe.keaki.s1@dc.tohoku.ac.jp", 
                                  notification:str = "abe"
                                  ):
        """
        MASAMUNEでjobを実行する際のスクリプトを作成する
        task...MDかDFT
        """
        def hours_to_hms(hours:float):
            total_seconds = int(hours * 3600)
            hh = total_seconds // 3600
            mm = (total_seconds % 3600) // 60
            ss = total_seconds % 60
            return f"{hh:02}:{mm:02}:{ss:02}"
        
        if walltime==None:
            if que=="P_016":
                walltime=72.0
            elif que=="DP_002":
                walltime=0.5
        
        lines:list = []
        # header
        lines += ["#!/bin/bash\n\n"]
        
        # options
        lines += [
            f"#PBS -l select={nodes}\n", 
            f"#PBS -l walltime={hours_to_hms(walltime)}\n",
            f"#PBS -q {que}\n",
            f"#PBS -N {jobname}_{task}\n",
            f"#PBS -M {emailaddress}\n",
            f"#PBS -m {notification}\n\n"
        ]
        
        # modules
        lines += [
            "# Cray環境からIntel環境に切り替え\n", 
            "module switch PrgEnv-cray PrgEnv-intel\n\n", 
            "# Crayの科学計算ライブラリをアンロード\n", 
            "module unload cray-libsci\n\n"
        ]
        
        # command
        if task.lower() == "md" or "molecular dynamics":
            software = "laich"
        elif task.lower() == "dft" or "density functional theory":
            software = "vasp_std"
            
        lines += [f"aprun -n {nodes*36} -N 36 -j 1 {software} > out.txt"]
        
        
        
##################################################
##################################################
##################################################
class ExportFile(ExportAtom, ExportPara, ExportMisc):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################