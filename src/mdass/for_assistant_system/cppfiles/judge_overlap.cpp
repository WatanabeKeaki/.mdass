#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // std::vectorの変換をサポートするために必要
#include <iostream>
#include <cassert> // assert マクロを使用するために必要
#include <math.h>

namespace py = pybind11;

float culc_dist(const std::vector<float> &point0, const std::vector<float> &point1)
{
    size_t numRows0 = point0.size();
    size_t numRows1 = point1.size();
    if (numRows0 != 3 || numRows1 != 3)
        abort();

    float temp = 0;
    for (int i = 0; i < 3; i++)
    {
        temp += pow(point0[i] - point1[i], 2);
    }
    float dist = pow(temp, 0.5);

    return dist;
}

std::vector<float> judge_overlap(const std::vector<std::vector<float>> &ixyz_remain,
                                 const std::vector<std::vector<float>> &ixyz_judge,
                                 const float &cutoff = 2.0)
{

    // 行列の行数と列数
    size_t numRows_remain = ixyz_remain.size();
    size_t numRows_judge = ixyz_judge.size();

    size_t numCols_remain = (numRows_remain > 0) ? ixyz_remain[0].size() : 0;
    size_t numCols_judge = (numRows_judge > 0) ? ixyz_judge[0].size() : 0;
    assert(numCols_remain == 4 || numCols_judge == 4);

    std::vector<float> selected_ids = {};

    // 行列の処理
    for (size_t i = 0; i < numRows_remain; ++i)
    {
        // 2から4番目の要素をfloatのベクトルに取り出す、xyz座標
        std::vector<float> xyz_remain(ixyz_remain[i].begin() + 1, ixyz_remain[i].begin() + 4);

        for (size_t j = 0; j < numRows_judge; ++j)
        {
            // 2から4番目の要素をfloatのベクトルに取り出す、xyz座標
            std::vector<float> xyz_judge(ixyz_judge[j].begin() + 1, ixyz_judge[j].begin() + 4);

            float dist = culc_dist(xyz_remain, xyz_judge);
            if (dist <= cutoff)
            {
                // 最初の要素をintに取り出す
                float id_judge = ixyz_judge[j][0];
                selected_ids.push_back(id_judge);
            }
        }
    }
    return selected_ids;
}

PYBIND11_MODULE(judge_overlap, m)
{
    m.def("judge_overlap", &judge_overlap, "judge atomic overlap");
}