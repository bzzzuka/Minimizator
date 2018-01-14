#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "sobolseqgenerator.h"

int SobolSeqGenerator::Init(uint32_t N_, uint32_t D_, std::string dir_file){
    N = N_;
    D = D_;
    current_point_number = -1;
    params_filename = dir_file;

    L = (unsigned)ceil(log((double)N)/log(2.0));

    std::ifstream infile(params_filename,std::ios::in);
    if (!infile) {
        std::cout << "Input file containing direction numbers cannot be found!" << std::endl;
        exit(1);
    }

    char buffer[1000];
    infile.getline(buffer,1000,'\n');

    dir_num_params.push_back(DirectionalNumbersParams());

    for (unsigned i = 1; i < D; i++)
    {
        dir_num_params.push_back(DirectionalNumbersParams());
        // Чтение параметров из файла
        unsigned d_, s_, a_;
        infile >> d_ >> s_ >> a_;
        dir_num_params[i].d = d_;
        dir_num_params[i].s = s_;
        dir_num_params[i].a = a_;
        dir_num_params[i].m.push_back(1);
        for (unsigned j = 1; j <= s_; j++)
        {
            unsigned m_;
            infile >> m_;
            dir_num_params[i].m.push_back(m_);
        }
    }
    infile.close();

    return 0;
}

PointReal SobolSeqGenerator::GeneratePoint(){
    int i,j,k;

    if (N == 0 || D == 0){
        std::cout << "Wrong initialization! N = " << N << ", D = " << D << std::endl;
        return PointReal();
    }

    if(current_point_number == N - 1){
        std::cout << "All " << N << " points are already generated!" << std::endl;
        return PointReal();
    }

    current_point_number++;

    if(current_point_number == 0)
    {
        //The 0-th point is zero
        last_generated_point = PointUnsigned(D, std::vector<unsigned>(D, 0));
        return PointReal(D, std::vector<Real>(D, 0));
    }

    PointUnsigned temp_point = PointUnsigned(D);

    C = 1;
    unsigned value = current_point_number - 1;
    while (value & 1) {
        value >>= 1;
        C++;
    }

    // Compute direction numbers V, scaled by pow(2,32)
    int v = 1 << (32 - C);
    temp_point.coordinate[0] = last_generated_point.coordinate[0] ^ v;

    // ----- Compute the remaining dimensions -----
    for (i = 1; i < D; i++){
        unsigned s = dir_num_params[i].s;
        unsigned a = dir_num_params[i].a;
        auto m     = dir_num_params[i].m;
        // Compute direction numbers V[1] to V[L], scaled by pow(2,32)
        std::vector<unsigned> V(L+1);

        if (L <= s){
            for (j = 1; j <= L; j++){
                V[j] = m[j] << (32 - j);
            }
        }
        else {
            for (j = 1; j <= s; j++){
                V[j] = m[j] << (32 - j);
            }
            for (j = s + 1; j <= L; j++) {
                V[j] = V[j - s] ^ (V[j - s] >> s);
                for (k = 1; k <= s - 1; k++)
                    V[j] ^= (((a >> (s - 1 - k)) & 1) * V[j - k]);
            }
        }

        temp_point.coordinate[i] = last_generated_point.coordinate[i] ^ V[C];
    }

    last_generated_point = temp_point;

    PointReal final_result = PointReal(D);
    for (i = 0; i < D; i++){
        final_result.coordinate[i] = (Real) temp_point.coordinate[i] / std::pow(Real(2.0), 32);
    }

    return final_result;
}