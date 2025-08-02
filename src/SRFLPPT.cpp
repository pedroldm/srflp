#include "SRFLPPT.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

SRFLPPT::SRFLPPT(std::string filename, int movementType,
                 double maxTempProportion)
    : mersenne_engine(rng_device()) {
    std::ifstream input(filename);
    if (!input.is_open()) {
        throw std::runtime_error("Error opening file: " + filename);
    }

    input >> this->n;

    this->lengths.resize(n);
    this->halfLengths.resize(n);
    this->frequencyMatrix.resize(n, std::vector<int>(n, 0));

    for (int i = 0; i < n; ++i) {
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
    SRFLPS ss;
    ss.cost = 0.0;
    ss.kl = {0, 0};
    ss.sol.resize(n);
    ss.ftfC.resize(n, std::vector<double>(n, 0.0));
    ss.ftfD.resize(n, std::vector<double>(n, 0.0));
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::shuffle(ss.sol.begin(), ss.sol.end(), this->mersenne_engine);
    return ss;
}

SRFLPS SRFLPPT::neighbor(SRFLPS sol) {
    SRFLPS newS = sol;
    std::uniform_int_distribution<> dist(0, this->n - 1);
    int index = dist(this->mersenne_engine);
    int newIndex = dist(this->mersenne_engine);

    switch (this->movementType) {
        case 1: /* swap */
            std::swap(newS.sol[index], newS.sol[newIndex]);
            break;
        case 2: /* 2-opt */
            if (index > newIndex) std::swap(index, newIndex);
            std::reverse(newS.sol.begin() + index,
                         newS.sol.begin() + newIndex + 1);
            break;
        case 3: /* re-insertion */
            if (index < newIndex) newIndex -= 1;
            int element = newS.sol[index];
            newS.sol.erase(newS.sol.begin() + index);
            newS.sol.insert(newS.sol.begin() + newIndex, element);
            newS.kl = {index, newIndex};
    }

    return newS;
}

double SRFLPPT::evaluate(SRFLPS &s) {
    !s.cost ? this->completeEvaluation(s) : this->deltaEvaluation(s);
    return s.cost;
}

void SRFLPPT::completeEvaluation(SRFLPS &s) {
    double dist;

    for (int i = 1; i < this->n; i++) {
        dist = this->halfLengths[s.sol[i]] + this->halfLengths[s.sol[i - 1]];
        s.ftfC[s.sol[i]][s.sol[i - 1]] =
            dist * this->frequencyMatrix[s.sol[i]][s.sol[i - 1]];
        s.ftfC[s.sol[i - 1]][s.sol[i]] =
            dist * this->frequencyMatrix[s.sol[i - 1]][s.sol[i]];
        s.ftfD[s.sol[i - 1]][s.sol[i]] = dist;
        s.ftfD[s.sol[i]][s.sol[i - 1]] = dist;
    }

    for (int leap = 2; leap < this->n; leap++) {
        for (int i = leap; i < this->n; i++) {
            dist =
                s.ftfD[s.sol[i - leap]][s.sol[i - 1]] +
                (this->halfLengths[s.sol[i - 1]] + this->halfLengths[s.sol[i]]);
            s.ftfC[s.sol[i]][s.sol[i - leap]] =
                dist * this->frequencyMatrix[s.sol[i]][s.sol[i - leap]];
            s.ftfC[s.sol[i - leap]][s.sol[i]] =
                dist * this->frequencyMatrix[s.sol[i - leap]][s.sol[i]];
            s.ftfD[s.sol[i - leap]][s.sol[i]] = dist;
            s.ftfD[s.sol[i]][s.sol[i - leap]] = dist;
        }
    }

    s.cost = 0.0;
    for (const auto &row : s.ftfC) {
        for (double val : row) {
            s.cost += val;
        }
    }
}

void SRFLPPT::deltaEvaluation(SRFLPS &s) {}
