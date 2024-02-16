#ifndef HEADERFILE_ARSENALHEADER
#define HEADERFILE_ARSENALHEADER

int GetInputs(int argc, char const *argv[]);
int LoadEoS();
int GenerateGrids();
int GenerateHSeq();
int CalcPsin();
int ExportPsin();

#endif