#ifndef HEADERFILE_ARSENALHEADER
#define HEADERFILE_ARSENALHEADER

int GetInputs(int argc, char const *argv[]);
int LoadEoS();
int LoadPhiGaussPts();
int GenerateGrids();
int GenerateHSeq();
int CalcjTn();
int ExportjTn();

#endif