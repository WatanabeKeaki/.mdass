mass_dict:dict = {
    "H" :   1.00794,        "He":   4.002602,       "Li":   6.941,
    "Be":   9.012182,       "B" :   10.811,         "C" :   12.0107,
    "N" :   14.0067,        "O" :   15.9994,        "F" :   18.9984032,
    "Ne":   20.1797,        "Na":   22.98976928,    "Mg":   24.305,
    "Al":   26.9815386,     "Si":   28.0855,        "P" :   30.973762,
    "S" :   32.065,         "Cl":   35.453,         "Ar":   39.948,
    "K" :   39.0983,        "Ca":   40.078,         "Sc":   44.955912,
    "Ti":   47.867,         "V" :   50.9415,        "Cr":   51.9961,
    "Mn":   54.938045,      "Fe":   55.845,         "Co":   58.933195,
    "Ni":   58.6934,        "Cu":   63.546,         "Zn":   65.38,
    "Ga":   69.723,         "Ge":   72.64,          "As":   74.9216,
    "Se":   78.96,          "Br":   79.904,         "Kr":   83.798,
    "Rb":   85.4678,        "Sr":   87.62,          "Y" :   88.90585,
    "Zr":   91.224,         "Nb":   92.90638,       "Mo":   95.96,
    "Tc":   99,             "Ru":   101.07,         "Rh":   102.9055,
    "Pd":   106.42,         "Ag":   107.8682,       "Cd":   112.411,
    "In":   114.818,        "Sn":   118.71,         "Sb":   121.76,
    "Te":   127.6,          "I" :   126.90447,      "Xe":   131.293,
    "Cs":   132.9054519,    "Ba":   137.327,        "La":   138.90547,
    "Ce":   140.116,        "Pr":   140.90765,      "Nd":   144.242,
    "Pm":   145,            "Sm":   150.36,         "Eu":   151.964,
    "Gd":   157.25,         "Tb":   158.92535,      "Dy":   162.5,
    "Ho":   164.93032,      "Er":   167.259,        "Tm":   168.93421,
    "Yb":   173.054,        "Lu":   174.9668,       "Hf":   178.49,
    "Ta":   180.94788,      "W" :   183.84,         "Re":   186.207,
    "Os":   190.23,         "Ir":   192.217,        "Pt":   195.084,
    "Au":   196.966569,     "Hg":   200.59,         "Tl":   204.3833,
    "Pb":   207.2,          "Bi":   208.9804,       "Po":   210,
    "At":   210,            "Rn":   222,            "Fr":   223,
    "Ra":   226,            "Ac":   227,            "Th":   232.03806,
    "Pa":   231.03588,      "U" :   238.02891}

# para.rdの各パラメータ
colm_list:list = [
#atom
["r_sigma",  "Val",     "mass",     "r_vdw",        "D",        "Y",        "r_pi",    "Val_e", 
 "arfa",     "Y_omega", "Val_boc",  "p_ovun5",      "unused0",  "X",        "eta",     "hflag", 
 "r_pipi",   "p_lp2",   "unused1",  "p_boc4",       "p_boc3",   "p_boc5",   "unused2", "unused3", 
 "p_ovun2",  "p_val3",  "unused4",  "Val_angle",    "p_val5",   "unused5",  "unused6", "unused7" ],
#bond
["D_sigma",    "D_pi",    "D_pipi", "p_be1",    "p_bo5",    "v13cor",   "p_bo6",    "p_ovun1", 
 "p_be2",      "p_bo3",   "p_bo4",  "unused0",  "p_bo1",    "p_bo2",    "ovc",      "unused1" ],
#offdiagonal
["D",   "r_vdw",      "arfa",   "r_sigma",   "r_pi",    "r_pipi"                             ],
#angle
["theta",  "p_val1",    "p_val2",    "p_coa1", "p_val7",    "p_pen1",  "p_val4"              ],
#torsion
["V_1",     "V_2",       "V_3",    "p_tor1", "p_cot1",   "unused0", "unused1"                ],
#hydrogenbond
["r_hb0",   "r_hb1",     "r_hb2",     "r_hb3"                                                ]]


keyword_laich_config = (
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

# keyword_vasp_INCAR = (
#     ("ALGO", "GGA", "PREC")
#     ("ISMEAR", "SIGMA")
    
    
# )

keyword_vasp_INCAR = {
    # 基本的な設定
    "ALGO": 
        """
        default: Normal
        DFT計算のアルゴリズムを指定する。
        久保研究室ではNormalやFastをよく使う。
        """, 
    "GGA": 
        """
        default: POTCARファイル依存
        LDAまたはGGAなど交換相関エネルギー汎関数を指定する。
        久保研究室ではPEをよく使う。
        """, 
    "PREC":
        """
        default: Normal
        計算の"精度"を指定する。
        具体的には、
        ENCUT: 平面基底波のエネルギーカットオフ
        NGX,NGY,NGZ: 高速フーリエ変換のグリッド数
        NGXF,NGXY,NGXZ: 高速フーリエ変換の細かいグリッド数
        PORT: プロジェクターが実空間でどれだけ正確に表現されるかを決定
        を自動で決定してくれる。        
        推奨される設定はNormalかAccurate
        """,
    "KSPACE":
        """
        default: 0.5
        KPOINTSファイルがない時にk点の数を決定する。
        値が小さいほどk点の数は増える。
        単位は/Å
        """,
    "IBRION":
        """
        defaut: -1
        原子、イオンの更新・移動方法を決定する。
        0: 分子動力学が実行される。
        1: RMM-DIISによる構造最適化。高速だが初期配置に大きく依存する。
        2: 共役勾配アルゴリズムによる最適化。RMM-DIISに比べ収束しやすい。
        3: 減衰分子動力学による緩和。POTIMの設定が重要となる。
        通常は2が勧められる。3はNEBで用いられた。
        """,
    "ISIF":
        """
        default: 0 (IBRION=0, 分子動力学)
                 2 (上記以外)
        原子座標・セルの自由度を指定。
        応力テンソルの計算、原子・イオンの位置やセル形状・サイズの扱いを指定。
        設定値: 応力テンソル計算, 変化可能な要素
        0: なし, 原子位置
        1: トレースのみ, 原子位置
        2: あり, 原子位置(分子動力学・NVTで利用)
        3: あり, 原子位置・セル形状・セルサイズ(分子動力学・NPTで利用)
        4: あり, 原子位置・セル形状
        5: あり, 原子位置・セルサイズ
        """,
    "POTIM":
        """
        default: None  (IBRION=0)
                 0.5   (IBRION=1,2,3)
                 0.015 (IBRION=5,6)
        分子動力学のタイムステップか、構造最適化のステップ幅を指定。
        分子動力学(IBRION=0)の場合、これを指定しないと計算はクラッシュする。
        RMM-DIIS(IBTION=1)の場合、POTIMの値が重要となる。
        vasp5.1以降のバージョンでは、POTIMが不当に大きい値なら、0.015に強制的に設定される。
        """, 
    "EDIFF":
        """
        default: 1.00e-4
        電子 SC ループのグローバル ブレーク条件を指定。
        EDIFFはeV 単位で指定。
        正確な結果を得たいなら、1.00e-6が強く推奨される。
        場合によっては1.00e-7を使う場合もあるが、収束しづらい。
        IBRIONの値に関わらず(分子動力学でも構造最適化でも)重要となる。
        """,
    "NSW":
        """
        default: 0
        計算のトータルステップを指定する。
        分子動力学(IBRION=0)では計算のステップ数を指定する。これを指定しないとクラッシュする。
        構造最適化(IBRION!=0)では計算の最大ステップ数を指定する。
        久保研究室では200に設定される場合が多い。
        """, 
    "NELM":
        """
        default: 60
        ESCが収束するまでのループを何回まで許容するかを設定する。
        vasp wikiには60で十分とあるが、初期構造が雑だったりすると60では全く収束しない。
        久保研究室では200に設定される場合が多い。
        """, 
    # 出力ファイルの設定
    "LWAVE":
        """
        default: .NOT. LH5
        計算の最後にWAVECARを出力するかを決定する。
        .FALSE.なら出力されない。
        """, 
    "LCHARGE":
        """
        default: .NOT. LH5 
        計算時に電荷密度が記載されているCHGCAR, CHGが書き込まれるかを決定する。
        .FALSE.なら書き込まれない。
        """, 
    # 発展的な設定
    "LREAL":
        """
        default: .FALSE.
        投影演算子を実空間で評価するか逆空間で評価するかを決定。
        原子数30以上では、Autoに設定することが推奨される。
        """, 
    "LASPH":
        """
        default: .FALSE.
        PAW球内の勾配に関する非球状の寄与を考慮するかを決定する。
        遷移金属酸化物を扱うにあたり、正確なエネルギーを計算するには
        .TRUE.にし、設定を有効にしなければならない。
        TiO2等を扱う際はこの設定をオンにする。
        """, 
    # 元素・モデル依存の設定
    "ISMEAR": # SIGMAの設定が必須
        """
        default: 1
        スメアリングを使用するかを設定する。
        スメアリングがオフ(ISMEAR=0)なら、電子はエネルギー準位を低い順に埋める。
        スメアリングがオン(ISMEAR>0)なら、エネルギー準位の高い値にもフェルミディラック統計に従い電子が入る。
        どの程度まで高い準位に電子を入れるかは、後述のSIGMAにより設定。
        ISMEARは-5～-1にも設定可能。
        金属なら1以上の値を入れ、半導体・絶縁体では0か-5を使用する。
        半導体・絶縁体で1以上の値を使用すると間違った結果が出る事がある。
        """, 
    "SIGMA": # ISMEARとセット
        """
        default: 0.2
        スメアリングの幅をeV単位で指定する。
        SIGMAが大きい程、高い準位まで電子が入る。
        どのような値を使うかは、先行研究やvasp wikiに添付されがちな例を参考にするとよい。
        """, 
    "ISPIN": 
        """
        default: 1
        スピン偏極の指定に利用。磁性のある金属(Fe等)で設定の必要がある。
        1: 非スピン偏極計算を実施
        2: スピン偏極計算を実施
        MAGMOMと合わせることで、共線磁性を扱える(本データセットには未実装)
        """, 
    "ISYM":
        """
        default: 1(VASPがUSPPで実行される場合)
                 3(LHFCALC=.TRUEの場合)
                 2(上記以外)
        計算における対称性の処理方法を指定する。
        久保研究室では対称性は利用しない事が多い。
        ISYM=-1, 0の場合、対称性の使用はオフになる。
        ISYM=0は分子動力学(IBRION=0)でのみ指定できる。
        構造最適化ではISYM=-1を指定。(分子動力学でも-1は指定可能)
        """, 
    # 使用するコア数の設定
    "NCORE":
        """
        default: 1
        個々の軌道で動作するコアの数を指定する。
        NCORE = 2 ~ ソケット・ノードあたりのコア数
        に設定することが推奨される。
        """, 
    "NPAR":
        """
        default: コア数
        並列処理されるバンドの数を決定する。
        defaultでは1つの軌道が1つのコアにより処理される。
        NPARが1の場合、NCOREがコア数に設定される。
        """,
    "KPAR":
        """
        default: 1
        並列処理されるk点の数を指定する。
        """,
    # 構造最適化
    "EDIFFG":
        """
        default: EDIFF*10
        構造最適化ループのブレーク条件を指定する。
        数値が小さいほど計算は収束しづらく、正確な結果となる。
        分子動力学(IBRION=0)には適用されない。
        負の値ならすべての力の合計がEDIFFGの絶対値より小さいときに緩和が停止。
        1.00e-5が推奨されるが、計算の種類によっては1.0e-2も使われる。
        NEB計算では-1.0e-2がサンプルとして用いられた。
        """,
    # 分子動力学, IBRION=0
    "TEBEG":
        """
        default: 0
        分子動力学(IBRION=0)における最初の温度を指定する。
        """,
    "TEEND":
        """
        default: TEBEGで指定した温度
        分子動力学(IBRION=0)における最終温度を指定する。
        """, 
    "MDALGO":
        """
        default: 0
        分子動力学のアルゴリズムを指定する。
        0: 標準, NVE(SMASS=-3に設定する)
        1: Andersen サーモスタット, NVT(ANDERSEN_PROBの指定が必要)
           LBLUEOUT =.TRUEと設定すれば、自由エネルギー勾配を計算できる。
        2: Nose-Hoover サーモスタット, NVT(SMASSの指定が必要)
           LBLUEOUT =.TRUEと設定すれば、自由エネルギー勾配を計算できる。
        3: Langevin サーモスタット, NVT(ISIF=2, LANGEVIN_GAMMAの指定が必要)
                                   NpT(ISIF=3, LANGEVIN_GAMMA・LANGEVIN_GAMMA_L・PMASSの指定が必要)
           NpTの場合, PSTRESSの指定により外部圧力も定義可能
        4: Nose Hoover chains サーモスタット, NVT(ISIF=2, NHC_NCHAINS, NHC_PERIODの指定が必要)
        """,
        "PSTRESS":
            """
            
            """,
        "PMASS":
            """
            
            """,
    # NEB
    "IOPT":
        """
        default: 0
        最適化のアルゴリズムを指定する。
        0: IBRIONに指定されたアルゴリズム
        1: Limited-memory Broyden-Fletcher-Goldfarb-Shanno
        2: Conjugate Gradient
        3: Quick-Min
        4: Steepest Descent
        7: Fast Inertial Relaxation Engine
        8: Machine learning (PyAMFF)
        サンプルでは3が指定された。
        """, 
    "ICHAIN":
        """
        default: 0
        実行するメソッドを指定。0ならNEBが実行される。
        """,
    "ICLIMB":
        """
        default: .TRUE.
        山登り法の様子を画像出力するかどうかを指定。
        """,
    "SPRING":
        """
        default: -5.0
        ばね定数。正の値なら通常のEB法になる。
        負の値ならNEBが実行される。
        default値からの変更はおすすめされない。
        """,
    "IMAGES":
        """
        default: None
        中間構造の数を指定。
        中間構造の数は奇数の方が良い。
        """
}