#include "arsenal.hpp"

int main(int argc, const char *argv[])
{
    if (GetInputs(argc, argv) != 1) return -1;
    if (LoadEoS()             != 1) return -2;
    if (LoadPhiGaussPts()     != 1) return -3;
    if (GenerateGrids()       != 1) return -4;
    if (GenerateHSeq()        != 1) return -5;
    if (CalcjTn()             != 1) return -6;
    if (ExportjTn()           != 1) return -7;

    return 0;
}