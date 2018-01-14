#pragma once

#include <string>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "tmsnet.hpp"

typedef Point<uint32_t> PointUnsigned;


struct DirectionalNumbersParams
{
    unsigned d, s, a;
    std::vector<unsigned> m;
};

class SobolSeqGenerator : public TMSNet
{
public:
    int Init(uint32_t N_, uint32_t D_, std::string dir_file);
    PointReal GeneratePoint();

private:
    uint32_t N;
    uint32_t D;
    uint32_t current_point_number;
    uint32_t L;
    uint32_t C;

    //Params
    std::string params_filename;
    std::vector<DirectionalNumbersParams> dir_num_params;

    PointUnsigned last_generated_point;
};

// #pragma once

// #include <vector>
// #include <string>
// #include "tmsnet.hpp"

// typedef Point<uint32_t> PointUnsigned;

// struct DirectionalNumbersParams
// {
//         unsigned d, s;
//         unsigned a;
//         std::vector<unsigned> m;
// };

// class SobolSeqGenerator : public TMSNet
// {
// public:
//     int Init();
//     int Init(uint32_t N_, uint32_t D_, std::string dir_file);
//     PointReal GeneratePoint();

// private:
//     uint32_t N;
//     uint32_t D;
//     uint32_t current_point_number;

//     // L = максимальное число нужных битов
//     uint32_t L;
//     uint32_t C;

//     // Параметры направляющих чисел
//     std::string params_filename;
//     std::vector<DirectionalNumbersParams> dir_num_params;

//     // Последняя сгенерированная одиночная точка
//     PointUnsigned last_generated_point;

// };
