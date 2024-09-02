import subprocess as sp
import pathlib as Path
import os
import shutil
from socket import gethostname
import subprocess
import time



class ExecMd:
    def __init__(self) -> None:
        pass
##################################################
    def exec_laich(self, 
                    workdir_laich:str="./", 
                    fname_out:str="./out.txt",
                    judge_done:bool = True, 
                    wait_process:bool = True):
        """_summary_
        laichをpythonから実行する。

        Args:
            path_wd_laich (str): MDを実行する場所。configやpara, inputがあるところ
            fname_out (str): 出力するファイル名
            judge_done (bool): out.txtがpath_wd_laichにあれば、実行を取りやめる
            wait_process (bool): Trueなら実行が終わるまでプログラムを止める

        Returns:
            _type_: _description_
        """
        #実行時に使うCPUコア数を決定する
        divx = int(self.config_laich[self.keyword_groups[4][0]])
        divy = int(self.config_laich[self.keyword_groups[4][1]])
        divz = int(self.config_laich[self.keyword_groups[4][2]])
        use_cpus = divx*divy*divz
        # laichを実行するディレクトリのパス
        path_wd_laich:Path.Path = self.path_wd.joinpath(workdir_laich)
        # laich実行時に出力するファイルのパス
        path_outfile:Path.Path = path_wd_laich.joinpath(fname_out)
        # outファイルが残っているなら、実行を取りやめて強制終了
        if judge_done and path_outfile.exists():
            print("out_file exist!")
            return 1
        # 実行
        hostname:str = gethostname()
        if "kbox" in hostname:
            cmd:str = f"mpiexec -np {use_cpus} laich"
        elif "super" in hostname or "mom" in hostname:
            cmd = f'aprun -n {use_cpus} -N 36 -j 1 laich'
        with open(path_outfile, "w") as f:
            self.process_laich = sp.Popen(
                                cmd,
                                # stderr = f,
                                stdout = f,
                                cwd = path_wd_laich,
                                shell= True
                                )
        if wait_process:
            self.process_laich.wait()
        return 0
##################################################
    def exec_vasp(self, num_process:int=1, num_nodes:int=1, vasp_command:str="vasp_std", fname_out:str="./out.txt"):
        hostname:str = gethostname()
        if "kbox" in hostname:
            cmd = f'mpiexec -np {num_process} {vasp_command} > {fname_out}'
        elif "super" in hostname or "mom" in hostname:
            cmd = f'aprun -n {num_process} -N {num_process / num_nodes} -j 1 {vasp_command} > {fname_out}'
        vasp_md_process = subprocess.Popen(cmd, cwd=self.path_wd, shell=True)
        time.sleep(5)
        return 0
##################################################
    def exec_packmol(self,
                     instance_bulk = None, 
                     instance_mole = None, 
                     pack_num:int = 1, 
                     tolerance:float = 2.0,
                     namedir_tmp:str = "packmol_tmp",
                     xyz_condition:list = None,
                     margin:float = 2.0, 
                     seed:int = -1,
                     packmol_cmd:str = f"packmol < {'packmol_mixture_comment.inp'}",
                     mask_mole:int = 0
                     ):
        pathdir_temp:Path.Path = self.path_wd.joinpath(namedir_tmp)
        os.makedirs(pathdir_temp, exist_ok=True)
        
        namefile_comment:str = "packmol_mixture_comment.inp"
        namefile_bulk:str = "bulk.xyz"
        namefile_mole:str = "mole.xyz"
        namefile_result:str = "packmol_mixture_result.xyz"
        
        pathfile_condition:Path.Path = pathdir_temp.joinpath(namefile_comment)
        pathfile_bulk:Path.Path = pathdir_temp.joinpath(namefile_bulk)
        pathfile_mole:Path.Path = pathdir_temp.joinpath(namefile_mole)
        pathfile_result:Path.Path = pathdir_temp.joinpath(namefile_result)
        
        if instance_bulk != None:
            atom_origin:dict = instance_bulk.atom.copy()
            num_atom_origin:int = instance_bulk.get_num_atom()
            cell_origin:dict = instance_bulk.cell.copy()
        
        if xyz_condition != None:
            cell_min:list = (cell_origin["min"]+margin).tolist()
            cell_max:list = (cell_origin["max"]-margin).tolist()
            xyz_condition:list = cell_min + cell_max
        else:
            assert len(xyz_condition) == 6
        
        # make condition_file
        lines:list = [
            f"tolerance {tolerance}\n",
            f"filetype xyz\n",
            f"seed {seed}\n",
            f"output packmol_mixture_result.xyz\n\n"
        ]
        
        if instance_bulk != None:
            lines += [
                "structure "+namefile_bulk+"\n",
                "\tnumber 1\n",
                "\tfixed 0. 0. 0. 0. 0. 0. 0. \n",
                "end structure\n\n",
            ]
        
        lines += [
            f"structure "+namefile_mole+"\n",
            f"\tnumber {pack_num}\n",
            f"\tinside box\t"+"\t".join(map(str, xyz_condition))+"\n",
            "end structure\n\n"
        ]
        
        with open(pathfile_condition, "w") as f:
            f.writelines(lines)
        
        # bulkのxyzファイルを出力
        if instance_bulk != None:
            instance_bulk.export_xyz(fname=pathfile_bulk, packmol=True)
        # moleculeのxyzファイルを出力
        instance_mole.export_xyz(fname=pathfile_mole, packmol=True)
        
        # packmol実行
        process = subprocess.Popen(packmol_cmd, cwd=pathdir_temp, shell=True)
        print(process)
        process.wait()
        
        #結果をインポート
        self.cell = instance_bulk.cell.copy()
        self.import_xyz(fname=pathfile_result, 
                        columns=["type","x","y","z"], 
                        sep="\s+"
                        )
        print(self.atom["id"])
        
        if instance_bulk != None:
            # atom
            for colm in self.colm_atom:
                self.atom[colm][0:num_atom_origin] = atom_origin[colm]
            self.atom["mask"][num_atom_origin:] = int(mask_mole)
            self.reset_index()
            # cell
            self.cell = cell_origin.copy()
             # mass
            self.make_mass()
        
        return 0
##################################################
##################################################
class AutoExec:
    def __init__(self) -> None:
        pass
##################################################
    def auto_laich(self,
                   judge_done:bool = True,):
        """
        optionはinput.rd(原子位置を格納)の機能指定に利用(fixなど)\n
        config.rd, input.rd, para.rdが必要
        """
        namefiles_input = [
            "input.rd",
            "config.rd",
            "para.rd",
        ]
        namefiles_output = [
            "Energy.dat", 
            "press.dat", 
            "out.txt", 
            "bonds_ovito", 
        ]
        restart = "restart"
        dump = "dump."
        
        # 結果を格納するディレクトリを作る
        dir_result = "result"
        dir_dump = "dumps"
        pathdir_result = self.path_wd.joinpath(dir_result)
        pathdir_dump = pathdir_result.joinpath(dir_dump)
        os.makedirs(pathdir_dump, exist_ok=True)
        
        # result内のout.txtからlaich実行回数を示すindexを決定
        out = namefiles_output[2]
        outs = [str(out_.name) for out_ in list(pathdir_result.iterdir()) if out in str(out_)]
        if outs:
            existing_numbers = [int(file_out.split("_")[0]) for file_out in outs if file_out.split("_")[0].isdigit()]
            index = max(existing_numbers) + 1
            pass
        else:
            index = 0
            
        # input, config, paraを作成
        # self.export_rd(options = options)
        # self.export_config_laich()
        # self.export_para()
        for namefile_input in namefiles_input:
            shutil.copy2(self.path_wd.joinpath(namefile_input), pathdir_result.joinpath(f"{index:02}_{namefile_input}"))
            shutil.copy2(self.path_wd.joinpath(namefile_input), pathdir_dump.joinpath(namefile_input))
            pass
        
        # laich実行、laichの実行終了までプログラムは待機
        self.exec_laich(judge_done=judge_done, wait_process=True)
        
        # 出力ファイルの処理
        for namefile_output in namefiles_output:
            try:
                shutil.move(self.path_wd.joinpath(namefile_output), pathdir_result.joinpath(f"{index:02}_{namefile_output}"))
            except:
                pass
        pathfiles_dump:list = [path_file for path_file in list(self.path_wd.iterdir()) if dump in str(path_file)]
        for pathfile_dump in pathfiles_dump:
            shutil.move(pathfile_dump, pathdir_dump.joinpath(pathfile_dump.name))
        try:
            shutil.copytree(self.path_wd.joinpath(restart), pathdir_result.joinpath(f"{index:02}_{restart}"))
        except FileNotFoundError:
            pass
        return 0
##################################################
##################################################
##################################################
class Exec(ExecMd, AutoExec):
    def __init__(self) -> None:
        super().__init__()
        pass
##################################################
##################################################
##################################################