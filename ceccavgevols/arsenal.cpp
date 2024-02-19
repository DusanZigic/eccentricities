#include "arsenal.hpp"
#include "linearinterpolation.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <random>

static const double PI = 3.1415926535897932L;

static std::vector<int> nList;
static unsigned int m;
static size_t QMCPtsN;

int GetInputs(int argc, char const *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && ((inputs[0] == "-h") || (inputs[0] == "--h"))) {
		std::cout << "default values: --nList=1,2,3,4,5,6,7,8 --m=2 --QMCPtsN=5000000" << std::endl;
		return 0;
	}

    std::map<std::string, std::string> inputparams;
	for (auto in : inputs)
	{
 	   	std::string key = in.substr(0, in.find("="));
 	   	std::string::size_type n = 0; while ((n = key.find("-", n)) != std::string::npos) {key.replace(n, 1, ""); n += 0;} //replacing all '-'
		std::string val = in.substr(in.find("=")+1, in.length());
		inputparams[key] = val;
	}

    std::string nStr = "1,2,3,4,5,6,7,8"; if (inputparams.count("n")       > 0)    nStr =           inputparams["n"];
                   m =                 2; if (inputparams.count("m")       > 0)       m = std::stoi(inputparams["m"]);
             QMCPtsN =           5000000; if (inputparams.count("QMCPtsN") > 0) QMCPtsN = std::stoi(inputparams["QMCPtsN"]);
    
    for (const auto& c : nStr) {
        if (c != ',') {
            nList.push_back(int(c - '0'));
        }
    }

    return 1;
}

static std::vector<std::vector<double>> QMCPtsX(4), QMCPtsY(4);

static double HaltonSequence(long unsigned int index, int base)
{
	double f = 1.0;
	double res = 0.0;

	while (index > 0) {
		f = f / base;
		res += f * (index % base);
		index = index / base; // integer division
	}

	return res;
}

int GenerateHSeq()
{
	std::vector<int> base{11, 13, 17, 19, 23, 29, 31, 37};
	std::random_device rd;
	auto rng = std::default_random_engine{rd()};
	std::shuffle(base.begin(), base.end(), rng);

	for (size_t i=0; i<QMCPtsN; i++) {
		for (size_t j=0; j<4; j++) {
			QMCPtsX[j].push_back(HaltonSequence((i+1)*409, base[j]));
			QMCPtsY[j].push_back(HaltonSequence((i+1)*409, base[j+4]));
		}
	}

	return 1;
}

static int LoadEvol(int n, interpFun& edensint)
{
    std::string path_in = "avgrotevoln" + std::to_string(n) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open evolution input file for n=" + std::to_string(n) + ". Aborting..." << std::endl;
        return -1;
    }

    std::vector<double> tauPts, xPts, yPts, tempPts, edPts;

    std::string line; double buffer;

    while (std::getline(file_in, line))
    {
        if (line.at(0) == '#')
            continue;
        
        std::stringstream ss(line);
        ss >> buffer;  tauPts.push_back(buffer);
        ss >> buffer;    xPts.push_back(buffer);
        ss >> buffer;    yPts.push_back(buffer);
        ss >> buffer; tempPts.push_back(buffer);
        ss >> buffer;   edPts.push_back(buffer);
    }

    edensint.SetData(tauPts, xPts, yPts, edPts);

    file_in.close();

    return 1;
}

std::map<int, double> Psin, Epsn;

static int CalcPsin(int n, interpFun& edensint)
{
    double tau0 = edensint.domain()[0][0];
    std::vector<double> xDomain = edensint.domain()[1];
    std::vector<double> yDomain = edensint.domain()[2];

    double sinsum = 0.0L, cossum = 0.0L;

	double x, y, rho, phi;

    #pragma omp parallel for reduction(+:sinsum,cossum) private(x,y,rho,phi)
	for (size_t iq=0; iq<QMCPtsN; iq++)
	{
		x = xDomain[0] + (xDomain[1] - xDomain[0])*QMCPtsX[0][iq];
		y = yDomain[0] + (yDomain[1] - yDomain[0])*QMCPtsY[0][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		sinsum += std::pow(rho, static_cast<double>(m))*std::sin(static_cast<double>(n)*phi)*edensint.interp(tau0, x, y);

		x = xDomain[0] + (xDomain[1] - xDomain[0])*QMCPtsX[1][iq];
		y = yDomain[0] + (yDomain[1] - yDomain[0])*QMCPtsY[1][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		cossum += std::pow(rho, static_cast<double>(m))*std::cos(static_cast<double>(n)*phi)*edensint.interp(tau0, x, y);
	}

	sinsum *= ((xDomain[1] - xDomain[0])*(yDomain[1] - yDomain[0])/static_cast<double>(QMCPtsN));
	cossum *= ((xDomain[1] - xDomain[0])*(yDomain[1] - yDomain[0])/static_cast<double>(QMCPtsN));

	Psin[n] = 1.0L/static_cast<double>(n) * std::atan2(sinsum, cossum) + PI/static_cast<double>(n);

	return 1;
}

static int CalcEpsn(int n, interpFun& edensint)
{
    double tau0 = edensint.domain()[0][0];
    std::vector<double> xDomain = edensint.domain()[1];
    std::vector<double> yDomain = edensint.domain()[2];

    double numerator = 0.0L, denominator = 0.0L;

	double x, y, rho, phi;

    #pragma omp parallel for reduction(+:numerator,denominator) private(x,y,rho,phi)
	for (size_t iq=0; iq<QMCPtsN; iq++)
	{
		x = xDomain[0] + (xDomain[1] - xDomain[0])*QMCPtsX[2][iq];
		y = yDomain[0] + (yDomain[1] - yDomain[0])*QMCPtsY[2][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		numerator += std::pow(rho, static_cast<double>(m))*std::cos(static_cast<double>(n)*(phi - Psin[n]))*edensint.interp(tau0, x, y);

		x = xDomain[0] + (xDomain[1] - xDomain[0])*QMCPtsX[3][iq];
		y = yDomain[0] + (yDomain[1] - yDomain[0])*QMCPtsY[3][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		denominator += std::pow(rho, static_cast<double>(m))*edensint.interp(tau0, x, y);
	}

	numerator   *= ((xDomain[1] - xDomain[0])*(yDomain[1] - yDomain[0])/static_cast<double>(QMCPtsN));
	denominator *= ((xDomain[1] - xDomain[0])*(yDomain[1] - yDomain[0])/static_cast<double>(QMCPtsN));

    Epsn[n] = -1.0L*numerator/denominator;

    return 1;
}

int CalcEcc()
{
    for (const auto& n : nList)
    {
        interpFun edensInt;
        if (LoadEvol(n, edensInt) != 1) return -1;
        if (CalcPsin(n, edensInt) != 1) return -2;
        if (CalcEpsn(n, edensInt) != 1) return -3;
    }

    return 1;
}

int ExportEcc()
{
    std::ofstream file_out("eccavgevols.dat", std::ios_base::out);
    if (!file_out.is_open()) {
        std::cerr << "Error: unable to open output file. Aborting..." << std::endl;
        return -1;
    }

    file_out << "#";
    for (const auto& n : nList)
        file_out << std::setw(11) << "Psi_" << n << " ";
    file_out << "\n" << " ";
    for (const auto& n : nList)
        file_out << std::scientific << std::setw(11) << std::setprecision(6) << Psin[n] << " ";
    file_out << "\n";
    
    file_out << "#";
    for (const auto& n : nList)
        file_out << std::setw(11) << "Eps_" << n << " ";
    file_out << "\n" << " ";
    for (const auto& n : nList)
        file_out << std::scientific << std::setw(11) << std::setprecision(6) << Epsn[n] << " ";

    file_out.close();

    return 1;
}