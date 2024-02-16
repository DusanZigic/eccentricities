#include "arsenal.hpp"

int main(int argc, const char *argv[])
{
    if (GetInputs(argc, argv) != 1) return -1;
    if (LoadEoS()             != 1) return -2;
    if (GenerateGrids()       != 1) return -3;
    if (GenerateHSeq()        != 1) return -4;
    if (CalcPsin()            != 1) return -5;
    if (ExportPsin()          != 1) return -6;

    return 0;
}