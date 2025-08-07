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
    std::vector<double> halfLengths;
    std::vector<std::vector<int>> frequencyMatrix;

    int movementType = 1;
    double maxTemp = 0.3;

    std::random_device rng_device;
    std::mt19937 mersenne_engine;

    SRFLPPT(std::string filename, int movementType, double maxTempProportion);
    SRFLPS construction();
    SRFLPS neighbor(SRFLPS sol);
    double evaluate(SRFLPS &sol);
    void init_E(std::vector<double>& E,int N,std::vector<int>& pk);
    void update_E(std::vector<double>& E,int k,int l,std::vector<int>& pk);
    double Fast_Inc_Eval(std::vector<double>& E,std::vector<int>& pk,int k,int l);
    double get_C(int k,std::vector<double>& E);
    void Insert_new(std::vector<int>& a,int k,int l);
    double completeEval(SRFLPS s);
};

#endif