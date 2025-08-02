#ifndef SRFLPPT_HPP
#define SRFLPPT_HPP

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "PTAPI/include/Problem.h"
#include "SRFLPS.hpp"

class SRFLPPT : public Problem<SRFLPS> {
   public:
    int n;
    std::vector<int> lengths;
    std::vector<int> halfLengths;
    std::vector<std::vector<int>> frequencyMatrix;

    int movementType = 1;
    double maxTemp = 0.3;

    std::random_device rng_device;
    std::mt19937 mersenne_engine;

    SRFLPPT(std::string filename, int movementType, double maxTempProportion);
    SRFLPS construction();
    SRFLPS neighbor(SRFLPS sol);
    double evaluate(SRFLPS sol);
    void oldCompleteEvaluation(SRFLPS s);
    void newCompleteEvaluation(SRFLPS s);
    void deltaEvaluation(SRFLPS s);
};

#endif