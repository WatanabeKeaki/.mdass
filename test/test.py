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

def all_zeros(lst):
    return all(x == 0 for x in lst)

# 使用例
lst1 = [0, 0, 0]
lst2 = [0, 1, 0]

print(all_zeros(lst1))  # 出力: True
print(all_zeros(lst2))  # 出力: False




# d = {0:1, 1:2}
# ds = {i*1000:d.copy() for i in range(3)}

# ds[1000][1] += 1
# ds = pd.DataFrame(ds).T
# print(ds)