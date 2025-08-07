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
    ss.sol.resize(n);
    ss.E.resize(n);
    // ss.ftfC.resize(n, std::vector<double>(n, 0.0));
    // ss.ftfD.resize(n, std::vector<double>(n, 0.0));
    std::iota(ss.sol.begin(), ss.sol.end(), 0);
    std::shuffle(ss.sol.begin(), ss.sol.end(), this->mersenne_engine);
    this->init_E(ss.E, this->n, ss.sol);
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

            double variation = this->Fast_Inc_Eval(newS.E, newS.sol, index, newIndex);
            this->update_E(newS.E, index, newIndex, newS.sol);
            newS.cost += variation;
    }

    return newS;
}

double SRFLPPT::evaluate(SRFLPS& s) {
    return s.cost;
}

double SRFLPPT::completeEval(SRFLPS s) {
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

void SRFLPPT::init_E(std::vector<double>& E, int N, std::vector<int>& pk) {
    for (int i = 0; i < N; i++) {
        E[i] = 0;
        for (int j = 0; j < N; j++) {
            if (i > j) {
                E[i] -= this->frequencyMatrix[pk[i]][pk[j]];
            } else if (i < j) {
                E[i] += this->frequencyMatrix[pk[i]][pk[j]];
            }
        }
    }
}

void SRFLPPT::update_E(std::vector<double>& E, int k, int l,
                       std::vector<int>& pk) {
    if (k < l) {
        double sum = 0, sum1 = 0, sum2 = 0;
        double temp_l = E[k];
        int limit = l - 1, i;
        for (i = k; i < limit; i += 2) {
            sum1 += (2 * this->frequencyMatrix[pk[k]][pk[i + 1]]);
            E[i] = E[i + 1] + 2 * this->frequencyMatrix[pk[k]][pk[i + 1]];

            sum2 += (2 * this->frequencyMatrix[pk[k]][pk[i + 2]]);
            E[i + 1] = E[i + 2] + 2 * this->frequencyMatrix[pk[k]][pk[i + 2]];
        }
        sum = sum1 + sum2;
        for (; i < l; i++) {
            sum += (2 * this->frequencyMatrix[pk[k]][pk[i + 1]]);
            E[i] = E[i + 1] + 2 * this->frequencyMatrix[pk[k]][pk[i + 1]];
        }
        E[l] = temp_l - sum;
    } else if (k > l) {
        double sum = 0, sum1 = 0, sum2 = 0;
        double temp_l = E[k];
        int limit = l + 2, i;
        for (i = k; i >= limit; i -= 2) {
            sum1 += 2 * this->frequencyMatrix[pk[k]][pk[i - 1]];
            E[i] = E[i - 1] - 2 * this->frequencyMatrix[pk[k]][pk[i - 1]];
            sum2 += 2 * this->frequencyMatrix[pk[k]][pk[i - 2]];
            E[i - 1] = E[i - 2] - 2 * this->frequencyMatrix[pk[k]][pk[i - 2]];
        }
        sum = sum1 + sum2;
        for (; i >= l + 1; i--) {
            sum += 2 * this->frequencyMatrix[pk[k]][pk[i - 1]];
            E[i] = E[i - 1] - 2 * this->frequencyMatrix[pk[k]][pk[i - 1]];
        }
        E[l] = temp_l + sum;
    }
}

double SRFLPPT::Fast_Inc_Eval(std::vector<double>& E, std::vector<int>& pk,
                              int k, int l) {
    double vary = 0;
    if (k < l) {
        if (l > k + 1) {
            double test = 0;
            test -= E[k];
            double e_sum = 0;
            for (int i = k + 1; i <= l - 1; i++) {
                e_sum += E[i];
                test += 2 * this->frequencyMatrix[pk[k]][pk[i]];
                vary += (test * (this->halfLengths[pk[i]] +
                                 this->halfLengths[pk[i + 1]]));
            }

            e_sum += E[l];
            if (k == 0) {
                vary +=
                    ((e_sum + test + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                     this->halfLengths[pk[k]]);
                vary += ((test + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                         this->halfLengths[pk[l]]);
                vary += ((this->get_C(k, E) * (this->halfLengths[pk[k]] -
                                               this->halfLengths[pk[k + 1]]) +
                          e_sum * this->halfLengths[pk[k]]));
            } else {
                vary +=
                    ((e_sum + test + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                     this->halfLengths[pk[k]]);
                vary += ((test + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                         this->halfLengths[pk[l]]);
                vary += ((e_sum + E[k]) * this->halfLengths[pk[k]]);
                vary -= (E[k] * this->halfLengths[pk[k + 1]]);
            }
        } else {
            if (k != 0) {
                vary -= ((E[k] + E[k + 1]) *
                         (this->halfLengths[pk[l]] - this->halfLengths[pk[k]]));
                vary +=
                    ((-E[k] + E[l] + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                     (this->halfLengths[pk[l]] + this->halfLengths[pk[k]]));
            } else {
                double temp_C = this->get_C(l, E);
                vary -= ((temp_C) *
                         (this->halfLengths[pk[l]] - this->halfLengths[pk[k]]));
                vary += ((2 * (E[l] + this->frequencyMatrix[pk[k]][pk[l]]) -
                          temp_C) *
                         (this->halfLengths[pk[l]] + this->halfLengths[pk[k]]));
            }
        }
    } else {
        if (k > l + 1) {
            double test = 0;
            test += E[k];
            double e_sum = 0;
            for (int i = k - 2; i >= l; i--) {
                e_sum += E[i];
                test += 2 * this->frequencyMatrix[pk[k]][pk[i + 1]];
                vary += ((test) * (this->halfLengths[pk[i]] +
                                   this->halfLengths[pk[i + 1]]));
            }
            e_sum += E[k - 1];
            if (l == 0) {
                vary +=
                    ((this->get_C(k - 1, E) * (this->halfLengths[pk[l]] -
                                               this->halfLengths[pk[k]]) +
                      (2 * this->frequencyMatrix[pk[k]][pk[l]] + test - e_sum) *
                          this->halfLengths[pk[l]]));
                vary += ((2 * this->frequencyMatrix[pk[k]][pk[l]] + test -
                          e_sum - E[k]) *
                         this->halfLengths[pk[k]]);
                vary += (E[k] * this->halfLengths[pk[k - 1]]);
            } else {
                vary += ((test + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                         this->halfLengths[pk[l]]);
                vary += ((2 * this->frequencyMatrix[pk[k]][pk[l]] + test -
                          e_sum - E[k]) *
                         this->halfLengths[pk[k]]);
                vary -= (e_sum * this->halfLengths[pk[k]]);
                vary += (E[k] * this->halfLengths[pk[k - 1]]);
            }
        } else {
            int temp = k;
            k = l;
            l = temp;
            if (k != 0) {
                vary -= ((E[k] + E[k + 1]) *
                         (this->halfLengths[pk[l]] - this->halfLengths[pk[k]]));
                vary +=
                    ((-E[k] + E[l] + 2 * this->frequencyMatrix[pk[k]][pk[l]]) *
                     (this->halfLengths[pk[l]] + this->halfLengths[pk[k]]));
            } else {
                double temp_C = this->get_C(l, E);
                vary -= ((temp_C) *
                         (this->halfLengths[pk[l]] - this->halfLengths[pk[k]]));
                vary += ((2 * (E[l] + this->frequencyMatrix[pk[k]][pk[l]]) -
                          temp_C) *
                         (this->halfLengths[pk[l]] + this->halfLengths[pk[k]]));
            }
        }
    }
    return vary;
}

double SRFLPPT::get_C(int k, std::vector<double>& E) {
    double C = E[0];
    for (int i = 1; i <= k; i++) {
        C += E[i];
    }
    return C;
}

void SRFLPPT::Insert_new(std::vector<int>& a, int k, int l) {
    int temp1 = a[k];
    if (k > l) {
        for (int i = k; i > l; i--) {
            a[i] = a[i - 1];
        }
    } else {
        for (int i = k; i < l; i++) {
            a[i] = a[i + 1];
        }
    }
    a[l] = temp1;
}