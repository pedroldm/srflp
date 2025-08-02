#ifndef SRFLP_HPP
#define SRFLP_HPP

#include "PTAPI/include/Problem.h"

struct SRFLPS : public solution {
    std::vector<int> sol;
    std::vector<std::vector<double>> ftfC;
    std::pair<int, int> kl;
    double cost;

    SRFLPS(int n);
};

#endif