#ifndef HEADERFILE_POLYINTHEADER
#define HEADERFILE_POLYINTHEADER

#include <vector>

double LinearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata);
double LinearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);
double CubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata);
double CubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit);

#endif