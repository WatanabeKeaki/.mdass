# import time
# from multiprocessing import Pool, Process

# def nijou(inputs):
#     x = inputs
#     print('input: %d' % x)
#     time.sleep(2)
#     retValue = x * x
#     return(retValue)

# if __name__ == "__main__":

#     # Pool()を定義
#     p = Pool()
#     results = [0 for _ in range(72)]
#     # プロセスを2つ非同期で実行
#     for i in range(72):
#         results[i] = p.apply_async(nijou, args=[i])

#     # 1秒間隔で終了チェックして終了したら結果を表示
#     for result in results:
#         if result.ready():
#             break
    
#     for result in results:
#         print(result.get())

#     p.close()

# def all_zeros(lst):
#     return all(x == 0 for x in lst)

# # 使用例
# lst1 = [0, 0, 0]
# lst2 = [0, 1, 0]

# print(all_zeros(lst1))  # 出力: True
# print(all_zeros(lst2))  # 出力: False




# d = {0:1, 1:2}
# ds = {i*1000:d.copy() for i in range(3)}

# ds[1000][1] += 1
# ds = pd.DataFrame(ds).T
# print(ds)

# def hours_to_hms(hours):
#     total_seconds = int(hours * 3600)
#     hh = total_seconds // 3600
#     mm = (total_seconds % 3600) // 60
#     ss = total_seconds % 60
#     return f"{hh:02}:{mm:02}:{ss:02}"

# # 入力が72時間の場合
# input_hours = 72
# hms_format = hours_to_hms(input_hours)
# print(hms_format)

# def ignore_comment(line):
#     config_dict = {}
#     line = line
#     # '#' でコメントを切り捨て、両端の空白を削除
#     line = line.split('#')[0]
#     if line:  # 空行でないか確認
#         # '=' でキーと値を分割
#         key_value = line.split('=')
#         if len(key_value) == 2:  # 有効なキー=値のペア
#             key = key_value[0].strip()
#             value = key_value[1].strip()
#             config_dict[key] = value
#     return config_dict

# line = "aaa = aa                  # a a a a                   "
# print(ignore_comment(line))
# i = 5

# print(list(range(-1,2)))

# from mdass import AssistantSystem
# import matplotlib.pyplot as plt
# import numpy as np

# a = AssistantSystem()
# x = np.arange(0, 10, 0.0005)
# y = np.sin(x*np.pi)

# print(x.shape[0])
# y = np.array(a.culc_moving_average(y, interbal=10000, ignore_part="as"))
# plt.plot(x,y)
# plt.savefig("test.png")

import numpy as np

A = 10
B = 10
C = 10

a = 60
b = 60
c = 60

cell0 = np.array([[A,0,0],[0,B,0],[0,0,C]])
