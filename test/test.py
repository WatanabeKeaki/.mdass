import matplotlib.pyplot as plt

# 例としていくつかのデータを生成します
x = range(10)
y1 = [i**2 for i in x]
y2 = [i**3 for i in x]

fig, ax1 = plt.subplots()

# 1番目の軸にプロットします
ax1.plot(x, y1, 'g-', label='Y1 data')  # label パラメータを追加して凡例用のラベルを指定します
ax1.set_xlabel('X data')
ax1.set_ylabel('Y1 data', color='g')

# 2番目の軸を追加します
ax2 = ax1.twinx()
ax2.plot(x, y2, 'b-', label='Y2 data')  # label パラメータを追加して凡例用のラベルを指定します
ax2.set_ylabel('Y2 data', color='b')

# 凡例を追加します
ax1.legend(loc='upper right')  # 凡例の位置を指定します

plt.savefig("test.png")
