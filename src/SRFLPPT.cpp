#include "SRFLPPT.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

SRFLPPT::SRFLPPT(std::string filename, int movementType, double maxTempProportion) : mersenne_engine(rng_device()) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    input >> this->n;

    this->lengths.resize(n);
    this->halfLengths.resize(n);
    this->frequencyMatrix.resize(n, std::vector<int>(n, 0));

    for(int i = 0 ; i < n ; ++i) {
        input >> this->lengths[i];
        this->halfLengths[i] = this->lengths[i] / 2.0;
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
    SRFLPS ss(this->n);
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::shuffle(ss.sol.begin(), ss.sol.end(), this->mersenne_engine);
    return ss;
}

SRFLPS SRFLPPT::neighbor(SRFLPS sol) {
    SRFLPS newS = sol;
    std::uniform_int_distribution<> dist(0, this->n - 1);
    int index = dist(this->mersenne_engine);
    int newIndex = dist(this->mersenne_engine);

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
    if (!s.cost) 
        this->newCompleteEvaluation(s);
    else
        this->deltaEvaluation(s);
    return s.cost;
}

void SRFLPPT::newCompleteEvaluation(SRFLPS s) {
    double dist;
    for(int i = 1 ; i < this->n ; i++) {
        dist = this->halfLengths[i] + this->halfLengths[i - 1];
        s.ftfC[i][i-1] = dist * this->frequencyMatrix[i][i-1];
        s.ftfC[i-1][i] = dist * this->frequencyMatrix[i - 1][i];
    }
}

void SRFLPPT::deltaEvaluation(SRFLPS s) {
    
}

void SRFLPPT::oldCompleteEvaluation(SRFLPS s) {
    s.cost = 0;
    for (int i = 0; i < this->n; i++) {
        for (int j = i + 1; j < this->n; j++) {
            int facility_i = s.sol[i];
            int facility_j = s.sol[j];

            double distance = this->halfLengths[facility_i] + this->halfLengths[facility_j];

            for (int k = i + 1; k < j; ++k) {
                int facility_k = s.sol[k];
                distance += this->lengths[facility_k];
            }

            double cost = this->frequencyMatrix[facility_i][facility_j];
            s.cost += cost * distance;
        }
    }
}
