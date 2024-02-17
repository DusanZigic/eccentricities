#ifndef HEADERFILE_ARSENALHEADER
#define HEADERFILE_ARSENALHEADER

int GetInputs(int argc, char const *argv[]);
int LoadEoS();
int LoadPsin();
int GetMinEvolLength();
int GenerateGrids();
int AvgRotatedEvols();

#endif