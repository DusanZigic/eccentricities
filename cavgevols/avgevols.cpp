#include "arsenal.hpp"

int main(int argc, const char *argv[])
{
    if (GetInputs(argc, argv) != 1) return -1;
    if (LoadEoS()             != 1) return -2;
    if (LoadPsin()            != 1) return -3;
    if (GetMinEvolLength()    != 1) return -4;
    if (GenerateGrids()       != 1) return -5;
    if (AvgRotatedEvols()     != 1) return -6;
}