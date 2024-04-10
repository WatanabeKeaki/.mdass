#include <iostream>
#include <vector>
#include <string>
#include <sstream>

int main()
{
    // 与えられたベクトル
    std::vector<std::string> vec = {"69", "0.104", "0.000", "0.000", "a"};

    // 1つ目の値をintに変換
    int first_value;
    std::istringstream(vec.front()) >> first_value;

    // 最後の値をfloatに変換
    float last_value;
    std::istringstream(vec.front()) >> last_value;
    // last_value += 1.0;

    // 結果の出力
    std::cout << "First value (int): " << first_value << std::endl;
    std::cout << "Last value (float): " << last_value << std::endl;

    return 0;
}
