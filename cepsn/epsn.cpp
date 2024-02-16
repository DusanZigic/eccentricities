#include "arsenal.hpp"

int main(int argc, char const *argv[])
{
    if (GetInputs(argc, argv) != 1) return -1;
    if (LoadEoS()             != 1) return -2;
    if (GenerateGrids()       != 1) return -3;
    if (GenerateHSeq()        != 1) return -4;
    if (LoadPsin()            != 1) return -5;
    if (CalcEpsn()            != 1) return -6;
    if (ExportEpsn()          != 1) return -7;

    return 0;
}