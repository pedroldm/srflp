#include "SRFLPS.hpp"

SRFLPS::SRFLPS(int n) {
    this->cost = 0.0;
    this->kl = {0, 0};
    this->sol.resize(n);
    this->ftfC.resize(n, std::vector<double>(n, 0.0));
}