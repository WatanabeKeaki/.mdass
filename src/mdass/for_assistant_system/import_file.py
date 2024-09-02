import numpy as np
import pathlib as Path
import pandas as pd
import re
from ase.build import molecule
# from ase.build import molecule

##################################################
##################################################
class ImportAtom:
    def __init__(self) -> None:
        self.colm_cell:list = ("min", "max")
        self.colm_mass:list = ("type", "mass")
        self.colm_atom:list = ("id", "type", "mask", 
                               "x", "y", "z", 
                               "vx", "vy", "vz", 
                               "fx", "fy", "fz", 
                               "q")
        self.colm_press = ("Step", 
                           "CellX", "CellY", "CellZ", 
                           "StressXX", "StressYY", "StressZZ", 
                           "StressXY", "StressYZ", "StressZX")
        
        self.keys_int:list = ("id", "type", "mask")
        
        self.path_wd: Path
        
        self.cell: dict
        self.mass: dict
        self.atom: dict
        pass
##################################################
    def import_rd(self, fname:str="input.rd") -> dict:
        """_summary_\n
        .rd形式のlaichで使われるinput fileを読み込む。
        Args:
            fname (str, optional): _description_. Defaults to "input.rd".
        """
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
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
            elif "#mass" in line:
                num_mass:int = int(line[1])
                startrow_mass:int = row+1
            elif "#atom" in line:
                num_atom:int = int(line[1])
                startrow_atom:int = row+1
        # cellsize
        self.cell:dict = dict()
        dates_cell:np.array = np.array([cellx, celly, cellz]).T # .Tは転置
        for key, date in zip(self.colm_cell, dates_cell):
            self.cell[key] = np.array(date).astype(np.float64)
        # mass
        self.mass:dict = dict()
        dates_mass:np.array = np.loadtxt(fname=str(path_input), 
                                         skiprows=startrow_mass, 
                                         max_rows=num_mass, 
                                         unpack=True)
        for key, date in zip(self.colm_mass, dates_mass):
            self.mass[key] = date
        # atom
        self.atom:dict = dict()
        dates_atom:np.array = np.loadtxt(fname=str(path_input), 
                                         skiprows=startrow_atom, 
                                         max_rows=num_atom, 
                                         unpack=True)
        for key, date in zip(self.colm_atom, dates_atom):
            self.atom[key] = date
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
    def import_pos(self, fname:str="dump.pos.0") -> dict:
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.split() for line in f.readlines(1000)]
        # 行ごとに解析
        for row, line in enumerate(lines):
            if not line:
                pass
            elif "NUMBER OF ATOMS" in " ".join(line):
                num_atom:int = int(lines[row+1][0])
            elif "BOX BOUNDS" in " ".join(line):
                cell_num = 3
                startrow_cell = row+1
            elif "ATOMS id type" in " ".join(line):
                startrow_atom = row+1
                colm_atom = line[2:]
                break
        # cellsize
        self.cell:dict = dict()
        dates_cell:np.array = np.loadtxt(fname=str(path_input), 
                                         skiprows=startrow_cell, 
                                         max_rows=cell_num, 
                                         unpack=True)
        for key, date in zip(self.colm_cell, dates_cell):
            self.cell[key] = date
        # atom
        self.atom:dict = dict()
        dates_atom:np.array = np.loadtxt(fname=str(path_input), 
                                         skiprows=startrow_atom, 
                                         max_rows=num_atom, 
                                         unpack=True)
        sorted_indices = np.argsort(dates_atom[0])
        dates_atom = dates_atom[:, sorted_indices]
        for key, date in zip(colm_atom, dates_atom):
            self.atom[key] = date
        # mass
        self.make_mass()
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
    def import_car(self, fname:str):
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.split() for line in f.readlines()]
        for i, line in enumerate(lines):
            if "PBC=OFF" == line[0]:
                pbc:bool = False
                startrow_atom = i+3
            if "PBC" == line[0]:
                pbc:bool = True
                # cell
                self.cell:dict = dict()
                self.cell[self.colm_cell[0]] = np.array([0., 0., 0.])
                self.cell[self.colm_cell[1]] = np.array(line[1:4], dtype=np.float64)
                startrow_atom = i+1
            if "end" == line[0]:
                endrow_atom = i-1
                num_atom = endrow_atom-startrow_atom+1
                break
        #atom
        self.atom:dict = dict()
        dates_atom:pd.DataFrame = pd.read_csv(filepath_or_buffer=path_input,  sep="\s+", 
                                              skiprows=startrow_atom, nrows=num_atom, 
                                              header=None)
        dates_atom.columns = ["A", "x", "y", "z", "B", "C", "D", "element", "E"]
        #id
        self.atom[self.colm_atom[0]] = np.array([i+1 for i in range(num_atom)])
        #type
        elements_temp:list = dates_atom["element"].to_list()
        types_temp:np.array = np.array([self.elm_type[element] for element in elements_temp])
        self.atom[self.colm_atom[1]] = types_temp
        #mask
        self.atom[self.colm_atom[2]] = np.array([0]*num_atom)
        #xyz
        for colm in self.colm_atom[3:6]:
            self.atom[colm] = dates_atom[colm].values
        #velocity, force, charge
        for colm in self.colm_atom[6:]:
            self.atom[colm] = np.array([0.]*num_atom)
        #cell
        if not pbc:
            self.cell:dict = dict()
            self.cell[self.colm_cell[0]] = np.array([0., 0., 0.])
            self.cell[self.colm_cell[1]] = np.array([0., 0., 0.])
            for i, q in enumerate(self.colm_atom[3:6]):
                q_min = np.min(self.atom[q])
                self.atom[q] -= q_min
                q_max = np.max(self.atom[q])
                self.cell[self.colm_cell[1]][i] = q_max+5
        # mass
        self.make_mass()
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
    def import_xyz(self, 
                   fname:str, 
                   columns:list=[], 
                   sep:str="\s+"
                   ):
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.rstrip() for line in f.readlines()]
            
        properties:tuple = ('id:I:1',     'species:S:1', 'pos:R:3',
                            'velo:R:3',   'mask:R:1',    'force:R:3',
                            'charge:R:1', 'velo_mag:R:1')
        
        for i, line in enumerate(lines):
            if i == 0:
                num_atom:int = int(line)
            if "Lattice" in line:
                self.cell:dict = dict()
                lattice = re.search(r'"([^"]*)"', line).group(1).split()
                lattice:list = list(map(float, lattice))
                
                cellx:float = lattice[0]
                celly:float = lattice[4]
                cellz:float = lattice[8]
                cell_max:list = [cellx, celly, cellz]
                self.cell[self.colm_cell[0]] = np.array([0., 0., 0.])
                self.cell[self.colm_cell[1]] = np.array(cell_max)
            
            if "Properties" in line and not columns:
                columns = []
                if 'id:I:1' in line:
                    columns.append(self.colm_atom[0])
                if 'species:S:1' in line:
                    columns.append(self.colm_atom[1])
                if 'pos:R:3' in line:
                    columns += self.colm_atom[3:6]
                if 'velo:R:3' in line:
                    columns += self.colm_atom[6:9]
                if 'mask:R:1' in line:
                    columns.append(self.colm_atom[2])
                if 'force:R:3' in line:
                    columns.append(self.colm_atom[9:12])
                if 'charge:R:1' in line:
                    columns.append(self.colm_atom[12])
                if 'velo_mag:R:1' in line:
                    columns.append("coming_soon(velo_mag)")
            if i==2:
                startrow_atom = i
                break
        #atom
        self.atom:dict = dict()
        dates_atom:pd.DataFrame = pd.read_csv(filepath_or_buffer=path_input,  sep=sep, 
                                              skiprows=startrow_atom, nrows=num_atom, 
                                              header=None)
        dates_atom.columns = columns
        #id
        if self.colm_atom[0] in columns:
            id_temp:np.array = dates_atom[self.colm_atom[0]].values
            self.atom[self.colm_atom[0]] = id_temp
        else:
            self.atom[self.colm_atom[0]] = np.arange(num_atom)+1
        #type
        if self.colm_atom[1] in columns:
            # typeの部分が数値に変換可能か(typeで書かれているか、elemetnで書かれているか)判定
            type_elem_ser = dates_atom[self.colm_atom[1]]
            type_or_elem:bool = type_elem_ser.apply(self.judge_numerical_conversion).all()
            if not type_or_elem: #Trueならtype, Falseならelemetnt
                elements_temp:list = type_elem_ser.to_list()
                types_temp:np.array = np.array([self.elm_type[element] for element in elements_temp])
            else: #typeの場合の処理
                types_temp:np.array = type_elem_ser.values
            self.atom[self.colm_atom[1]] = types_temp
        #mask
        if self.colm_atom[2] in columns:
            mask_temp:np.array = dates_atom[self.colm_atom[2]].values
            self.atom[self.colm_atom[2]] = mask_temp.astype(int)
        #xyz
        for colm in self.colm_atom[3:6]:
            self.atom[colm] = dates_atom[colm].values
        #velocity
        if self.colm_atom[6] in columns:
            for colm in self.colm_atom[6:9]:
                self.atom[colm] = dates_atom[colm].values
        #force
        if self.colm_atom[9] in columns:
            for colm in self.colm_atom[9:12]:
                self.atom[colm] = dates_atom[colm].values
        #charge
        if self.colm_atom[12] in columns:
            mask_temp:np.array = dates_atom[self.colm_atom[12]].values
            self.atom[self.colm_atom[12]] = mask_temp.astype(int)
        # mass
        self.make_mass()
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
    def import_molecule(self, name_mole:str, margin:float=1.5):
        self.atom:dict=dict()
        self.cell:dict=dict()
        
        atom_temp_ = molecule(name_mole)
        
        # type
        elements_temp_:list = atom_temp_.get_chemical_symbols()
        types_temp_:np.array = np.array([self.elm_type[element] for element in elements_temp_])
        self.atom[self.colm_atom[1]] = types_temp_
        
        # id
        num_atom:int = types_temp_.shape[0]
        ids_temp_:np.array = np.arange(num_atom)+1
        self.atom[self.colm_atom[0]] = ids_temp_
        
        # mask
        masks_temp_ = np.zeros(num_atom, dtype=int)
        self.atom[self.colm_atom[2]] = masks_temp_
        
        # position
        positions_temp_:np.array = atom_temp_.positions.T
        lens_min:list = [None]*3 # cellの決定用
        for ind, colm_q in enumerate(self.colm_atom[3:6]):
            self.atom[colm_q] = positions_temp_[ind]
            q_min:float = np.min(positions_temp_[ind])
            lens_min[ind] = q_min
        lens_max:list = [None]*3
        for ind, colm_q in enumerate(self.colm_atom[3:6]):
            self.atom[colm_q] -= lens_min[ind]
            q_max:float = np.max(self.atom[colm_q])
            lens_max[ind] = q_max
        
        # velocity, force, charge
        zeros_array = np.zeros(num_atom)
        for colm in self.colm_atom[6:]:
            self.atom[colm] = zeros_array.copy()
        
        # cell
        self.cell[self.colm_cell[0]] = np.zeros(3)
        self.cell[self.colm_cell[1]] = np.array(lens_max)+margin
        
        # mass
        self.make_mass()
        
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
    
    def import_XDATCAR(self, fname:str="XDATCAR", step:int=None):
        """
        step=Noneなら、最後のstepが選択される
        """
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.rstrip() for line in f.readlines()]
        
        steps:int = 0
        for i, line in enumerate(lines):
            if "Direct configuration" in line:
                steps += 1
            if i==1:
                scaling_factor = float(line)
            if i==2:
                len_x = float(line.split()[0])*scaling_factor
            if i==3:
                len_y = float(line.split()[1])*scaling_factor
            if i==4:
                len_z = float(line.split()[2])*scaling_factor
            if i==5:
                elements_temp__:list = line.split()
            if i==6:
                numatom_per_type:list = list(map(int,line.split()))
                num_atom:int = sum(numatom_per_type)
        if step == None:
            step:int = steps
        
        # cell
        self.cell[self.colm_cell[0]] = np.zeros(3)
        self.cell[self.colm_cell[1]] = np.array([len_x,len_y,len_z])
        # atom
        # id
        self.atom[self.colm_atom[0]] = np.arange(num_atom)+1
        # type
        elements_temp_:list = []
        for num, element in zip(numatom_per_type, elements_temp__):
            elements_temp_ += [element]*num
            types_temp_:np.array = np.array([self.elm_type[element] for element in elements_temp_])
            self.atom[self.colm_atom[1]] = types_temp_
        # pos
        skiprows:int = 7+(num_atom+1)*(step-1)+1
        positions = np.loadtxt(path_input, skiprows=skiprows, max_rows=num_atom).T
        xs_temp_:np.array = positions[0]*len_x
        ys_temp_:np.array = positions[1]*len_y
        zs_temp_:np.array = positions[2]*len_z
        self.atom[self.colm_atom[3]] = xs_temp_
        self.atom[self.colm_atom[4]] = ys_temp_
        self.atom[self.colm_atom[5]] = zs_temp_
        
        # mass
        self.make_mass()
        
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
    
    def import_POSCAR(self, fname:str="POSCAR"):
        # パスの作成
        path_input = self.path_wd.joinpath(fname)
        # ファイルの読み込み
        with open(path_input, "r") as f:
            lines = [line.rstrip() for line in f.readlines()]
        
        for i, line in enumerate(lines):
            if i==1:
                scaling_factor = float(line)
            if i==2:
                len_x = float(line.split()[0])*scaling_factor
            if i==3:
                len_y = float(line.split()[1])*scaling_factor
            if i==4:
                len_z = float(line.split()[2])*scaling_factor
            if i==5:
                elements_temp__:list = line.split()
            if i==6:
                numatom_per_type:list = list(map(int,line.split()))
                num_atom:int = sum(numatom_per_type)
            if i == 7:
                mode:str = line.split()[0]
                skiprows:int = i
                break
            
        # cell
        self.cell[self.colm_cell[0]] = np.zeros(3)
        self.cell[self.colm_cell[1]] = np.array([len_x,len_y,len_z])
        # atom
        # id
        self.atom[self.colm_atom[0]] = np.arange(num_atom)+1
        # type
        elements_temp_:list = []
        for num, element in zip(numatom_per_type, elements_temp__):
            elements_temp_ += [element]*num
            types_temp_:np.array = np.array([self.elm_type[element] for element in elements_temp_])
            self.atom[self.colm_atom[1]] = types_temp_
        # pos
        positions = np.loadtxt(path_input, skiprows=skiprows, max_rows=num_atom).T
        xs_temp_:np.array = positions[0]
        ys_temp_:np.array = positions[1]
        zs_temp_:np.array = positions[2]
        self.atom[self.colm_atom[3]] = xs_temp_
        self.atom[self.colm_atom[4]] = ys_temp_
        self.atom[self.colm_atom[5]] = zs_temp_
        
        # mass
        self.make_mass()
        
        # correction
        self.correct_dates()
        self.correct_dtype()
        return 0
##################################################
##################################################
##################################################
class ImportPara:
    def __init__(self) -> None:
        pass
##################################################
    def import_para(self, fname:str="para.rd") -> dict:
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
        path_input = self.path_wd.joinpath(fname)
        # para.rdを読み込む為に使用
        checkword_list = ['Nr of atom', 'Nr of bonds','Nr of off-diagonal', 
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
    def import_bond(self, fname:str="dump.bond.0"):
        path_input = self.path_wd.joinpath(fname)
        
        with open(path_input, "r") as f:
            current_row = 0
            while current_row < 5:
                line = f.readline()
                if len(line) == 1:
                    continue
                if 'NUMBER' in line:
                    num_atom = int(f.readline())
                    current_row += 1 # readline()を行った為
                    skiprows = current_row+1
                    break
                current_row += 1
        
        #律速
        df_bond:pd.DataFrame = pd.read_csv(path_input, skiprows=skiprows,  names=[0,1,2,3,4],
                                           sep=" ", header=None, engine="c")
        df_bond:pd.DataFrame = df_bond.replace('Atom', -1).astype("float32")
        array_bond:np.array = df_bond.values
        del df_bond
        
        # bondorderを抽出
        bondorders:np.array = array_bond[:, 1:5]
        bondids = array_bond[:, 0].astype(int)
        
        # 欠損値を含む行を選択するマスクを作成する,bool型をvalueに持つnp.arrayが作成される
        mask_nan = np.any(np.isnan(array_bond), axis=1)
        rownums_mainhead = np.where(mask_nan)[0]
        
        # 欠損値を含む行の2番目と3番目の列を選択し、整数型に変換する
        mainid_bondnum = array_bond[mask_nan, 1:3].astype(int).T
        mainids = mainid_bondnum[0]
        bondnums = mainid_bondnum[1]
        
        # 第二律速
        # 以下2つのarrayは一対一対応でなければならない
        self.bondorders_list = [None]*num_atom
        self.bondids_list = [None]*num_atom
        # print(num_atom)
        
        for rownum, mainid, bondnum in zip(rownums_mainhead, mainids, bondnums):
            mainidx = mainid-1 # 以降, main_idxの定義
            # print(mainidx)
            self.bondorders_list[mainidx] = bondorders[rownum+1: rownum+bondnum+1]
            self.bondids_list[mainidx] = bondids[rownum+1: rownum+bondnum+1]
        return 0
##################################################
    def import_bond2(self, fname:str="dump.bond.0"):
        path_input = self.path_wd.joinpath(fname)
        
        with open(path_input, "r") as f:
            current_row = 0
            while current_row < 5:
                line = f.readline()
                if len(line) == 1:
                    continue
                if 'NUMBER' in line:
                    num_atom = int(f.readline())
                    current_row += 1 # readline()を行った為
                    skiprows = current_row+1
                    break
                current_row += 1
        
        #律速
        df_bond:pd.DataFrame = pd.read_csv(path_input, skiprows=skiprows,  names=[0,1,2,3,4],
                                           sep=" ", header=None, engine="c")
        df_bond:pd.DataFrame = df_bond.replace('Atom', -1).astype("float32")
        array_bond:np.array = df_bond.values
        del df_bond
        
        # bondorderを抽出
        bondorders:np.array = array_bond[:, 1:5]
        bondids = array_bond[:, 0].astype(int)
        
        # 欠損値を含む行を選択するマスクを作成する,bool型をvalueに持つnp.arrayが作成される
        mask_nan = np.any(np.isnan(array_bond), axis=1)
        rownums_mainhead = np.where(mask_nan)[0]
        
        # 欠損値を含む行の2番目と3番目の列を選択し、整数型に変換する
        mainid_bondnum = array_bond[mask_nan, 1:3].astype(int).T
        mainids = mainid_bondnum[0]
        bondnums = mainid_bondnum[1]
        
        # 第二律速
        # 以下2つのarrayは一対一対応でなければならない
        self.bondorders_list = [None]*num_atom
        self.bondids_list = [None]*num_atom
        # print(num_atom)
        
        for rownum, mainid, bondnum in zip(rownums_mainhead, mainids, bondnums):
            mainidx = mainid-1 # 以降, main_idxの定義
            # print(mainidx)
            self.bondorders_list[mainidx] = bondorders[rownum+1: rownum+bondnum+1]
            self.bondids_list[mainidx] = bondids[rownum+1: rownum+bondnum+1]
        return 0
##################################################
##################################################
class ImportMisc:
    def __init__(self) -> None:
        pass
##################################################
    def import_config_laich(self, fname:str="config.rd") -> dict:
        path_input = self.path_wd.joinpath(fname)
        with open(path_input, "r") as f:
            lines = f.readlines()
        self.config_laich:dict = dict()
        for line in lines:
            if line.strip() != "" and "#" not in line:
                key = line.split()[0]
                value = line.split()[1]
                self.config_laich[key] = value
        return 0
##################################################
    def import_INCAR(self, fname:str="INCAR"):
        path_input = self.path_wd.joinpath(fname)
        with open(path_input, "r") as f:
            lines:list = f.readlines()
        self.incar:dict = dict()
        for line in lines:
            key, value = map(str.strip, line.split('='))
            self.incar[key] = value
        return 0
##################################################
    def import_press_laich(self, fname:str="press.dat"):
        """
        laichから出力されるpress.datを読み込む。
        """
        path_input = self.path_wd.joinpath(fname)
        with open(path_input, "r") as f:
            line:str = f.readline()
        if tuple(line.split()) == self.colm_press:
            skiprows = 1
        else:
            skiprows = 0
        press_temp = np.loadtxt(fname=path_input, 
                                skiprows=skiprows, 
                                unpack=True)
        if not self.press:
            for colm, values in zip(self.colm_press, press_temp):
                self.press[colm] = values
        else:
            for colm, values in zip(self.colm_press, press_temp):
                origin = self.press[colm]
                add = values
                self.press[colm] = np.concatenate([origin, add])
            sort_indices = np.argsort(self.press[self.colm_press[0]])
            for colm, value in self.press.items():
                self.press[colm] = value[sort_indices]
        return 0
##################################################
##################################################
##################################################
class ImportFile(ImportAtom, ImportPara, ImportBond, ImportMisc):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################