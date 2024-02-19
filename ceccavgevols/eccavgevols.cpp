#include "arsenal.hpp"

int main(int argc, char const *argv[])
{
    if (GetInputs(argc, argv) != 1) return -1;
    if (GenerateHSeq()        != 1) return -2;
    if (CalcEcc()             != 1) return -3;
    if (ExportEcc()           != 1) return -4;

    return 0;
}