#include "SRFLPPT.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

SRFLPPT::SRFLPPT(std::string filename, int movementType, double maxTempProportion) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    input >> this->n;

    this->lengths.resize(n);
    this->frequencyMatrix.resize(n, std::vector<int>(n, 0));

    for(int i = 0 ; i < n ; ++i) {
        input >> this->lengths[i];
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            input >> this->frequencyMatrix[i][j];
        }
    }

    this->movementType = movementType;
    this->maxTemp = this->n * maxTempProportion;

    input.close();
}

SRFLPS SRFLPPT::construction() {
    SRFLPS ss;
    ss.sol.resize(this->n);
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::random_device rnd_device;
    std::mt19937 mersenne_engine{rnd_device()};
    std::shuffle(ss.sol.begin(), ss.sol.end(), mersenne_engine);
    return ss;
}

SRFLPS SRFLPPT::neighbor(SRFLPS sol) {
    SRFLPS newS = sol;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dist(0, this->n - 1);
    int index = dist(gen);
    int newIndex = dist(gen);

    switch(this->movementType) {
        case 1 : /* swap */
            std::swap(newS.sol[index], newS.sol[newIndex]);
            break;
        case 2 : /* 2-opt */
            if (index > newIndex) 
                std::swap(index, newIndex);
            std::reverse(newS.sol.begin() + index, newS.sol.begin() + newIndex + 1);
            break;
        case 3 : /* re-insertion */
            if (index < newIndex)
                newIndex -= 1;
            int element = newS.sol[index];
            newS.sol.erase(newS.sol.begin() + index);
            newS.sol.insert(newS.sol.begin() + newIndex, element);
    }

    return newS;
}

double SRFLPPT::evaluate(SRFLPS s) {
    double totalCost = 0;
    for (int i = 0; i < this->n; i++) {
        for (int j = i + 1; j < this->n; j++) {
            int facility_i = s.sol[i];
            int facility_j = s.sol[j];

            double distance = this->lengths[facility_i] / 2.0 + this->lengths[facility_j] / 2.0;

            for (int k = i + 1; k < j; ++k) {
                int facility_k = s.sol[k];
                distance += this->lengths[facility_k];
            }

            double cost = this->frequencyMatrix[facility_i][facility_j];
            totalCost += cost * distance;
        }
    }
    return totalCost;
}
