#include "PTAPI/include/PT.h"
#include "SRFLPPT.hpp"

#include <iostream>
#include <string>
#include <sstream>

int main(int argc, char* argv[]) {
    float tempMin = 0.3f;
    float maxTempProportion = 1.0f;
    int tempL = 12;
    float MKL = 400;
    int PTL = 2000;
    int tempD = 4;
    int upType = 1;
    int tempUpdate = 3;
    int movementType = 3;
    string filePath;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.find("--tempMin=") == 0) {
            std::istringstream(arg.substr(10)) >> tempMin;
        } else if (arg.find("--tempL=") == 0) {
            std::istringstream(arg.substr(8)) >> tempL;
        } else if (arg.find("--MKL=") == 0) {
            std::istringstream(arg.substr(6)) >> MKL;
        } else if (arg.find("--PTL=") == 0) {
            std::istringstream(arg.substr(6)) >> PTL;
        } else if (arg.find("--tempD=") == 0) {
            std::istringstream(arg.substr(8)) >> tempD;
        } else if (arg.find("--upType=") == 0) {
            std::istringstream(arg.substr(9)) >> upType;
        } else if (arg.find("--tempUpdate=") == 0) {
            std::istringstream(arg.substr(13)) >> tempUpdate;
        } else if (arg.find("--maxTempProportion=") == 0) {
            std::istringstream(arg.substr(20)) >> maxTempProportion;
        } else if (arg.find("--movementType=") == 0) {
            std::istringstream(arg.substr(15)) >> movementType;
        } else if (arg.find("--filePath=") == 0) {
            std::istringstream(arg.substr(11)) >> filePath;
        } else {
            throw std::runtime_error("Unkown argument: " + arg);
        }
    }

    SRFLPPT* prob = new SRFLPPT(filePath, movementType, maxTempProportion);
    PT<SRFLPS> algo(
        tempMin,
        prob->maxTemp,
        tempL,
        MKL,
        PTL,
        tempD,
        upType,
        std::max(PTL / tempUpdate, 1)
    );
    SRFLPS sol = algo.start(15, prob);
    std::cout << std::fixed << std::setprecision(2) << sol.evalSol << std::endl;    
    delete prob;
    return 0;
}