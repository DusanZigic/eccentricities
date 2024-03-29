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

static unsigned int m; static size_t QMCPtsN;

int GetInputs(int argc, char const *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && ((inputs[0] == "-h") || (inputs[0] == "--h"))) {
		std::cout << "default values: --eventN=1000 --m=2 --QMCPtsN=5000000" << std::endl;
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

          m =       2; if (inputparams.count("m")       > 0)       m = std::stoi(inputparams["m"]);
    QMCPtsN = 5000000; if (inputparams.count("QMCPtsN") > 0) QMCPtsN = std::stoi(inputparams["QMCPtsN"]);

    return 1;
}

static interpFun EoS;

int LoadEoS()
{
	std::string path_in = "./eos/eos.dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open EoS file. Aborting..." << std::endl;
		return -1;
	}

    std::vector<double> temps; temps.push_back(0.0);
	std::vector<double> edens; edens.push_back(0.0);

	std::string line; double buff_e, buff_p, buff_s, buff_T;

	while (std::getline(file_in, line))
	{
		if (line.rfind("#", 0) == 0) continue;

		std::stringstream ss(line);
		ss >> buff_e;
		ss >> buff_p;
		ss >> buff_s;
		ss >> buff_T;

		temps.push_back(buff_T); edens.push_back(buff_e);
	}

	file_in.close();

	EoS.SetData(temps, edens);

	return 1;
}

static std::vector<double> xGridPts, yGridPts;

int GenerateGrids()
{
	std::string path_in = "./tempevols/temp_grids.dat";
	std::ifstream file_in(path_in, std::ios_base::in);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open grid parameters file. Aborting..." << std::endl;
		return -1;
	}

	double tau_min, tau_max, tau_step; int tau_step_n = 0;
	double   x_min,   x_max,   x_step; int   x_step_n = 0;
	double   y_min,   y_max,   y_step; int   y_step_n = 0;

	std::string line; double buff;

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buff; tau_min  = buff;
							   		ss >> buff; tau_step = buff;
							   			        tau_max  = 25.0;
		double tau = tau_min; while (tau < (tau_max+tau_step)) {tau_step_n++; tau+=tau_step;}
	}

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buff; x_min  = buff;
							   		ss >> buff; x_max  = buff;
							   		ss >> buff; x_step = buff;
		double x = x_min; while (x < (x_max+x_step)) {x_step_n++; x+=x_step;}
	}

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buff; y_min  = buff;
							   		ss >> buff; y_max  = buff;
							   		ss >> buff; y_step = buff;
		double y = y_min; while (y < (y_max+y_step)) {y_step_n++; y+=y_step;}
	}

	file_in.close();

	for (int i_x=0; i_x<x_step_n; i_x++) xGridPts.push_back(x_min + i_x*x_step);
	for (int i_y=0; i_y<y_step_n; i_y++) yGridPts.push_back(y_min + i_y*y_step);

	return 1;
}

static std::vector<std::vector<double>> QMCPtsX(2), QMCPtsY(2);

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
		for (size_t j=0; j<2; j++) {
			QMCPtsX[j].push_back(HaltonSequence((i+1)*409, base[j]));
			QMCPtsY[j].push_back(HaltonSequence((i+1)*409, base[j+2]));
		}
	}

	return 1;
}

static std::map<size_t, std::map<int, double>> Psin;
static size_t eventN = 0;
static std::vector<int> nList;

int LoadPsin()
{
    std::string path_in = "psin_m=" + std::to_string(m) + ".dat";;
    std::ifstream file_in(path_in, std::ios_base::in);
    if(!file_in.is_open()) {
        std::cerr << "Error: unable to open psin file. Aborting..." << std::endl;
        return -1;
    }

    std::string line; size_t event_id; std::string sbuffer; double dbuffer;

	while (std::getline(file_in, line))
	{
		if (line.at(0) == '#') {
			std::stringstream ss(line);
            while (ss >> sbuffer) {
                if (sbuffer.find("Psi") != std::string::npos)
                    nList.push_back(std::stoi(sbuffer.substr(sbuffer.find("_") + 1, sbuffer.length())));
            }
			continue;
		}
		
		std::stringstream ss(line);
		ss >> event_id;
		for (const auto& n : nList) {
			ss >> dbuffer;
			Psin[event_id][n] = dbuffer;
		}
	}

    file_in.close();

    eventN = Psin.size();

    return 1;
}

static int LoadEvol(size_t event_id, interpFun& edensint)
{
	std::string path_in = "tempevols/TProfile_" + std::to_string(event_id) + ".dat";
	std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
	if (!file_in.is_open()) {
		std::cerr << "Error: unable to open evolution file for event " << event_id << ". Aborting..." << std::endl;
	}

	std::vector<std::vector<double>> edens(xGridPts.size(), std::vector<double>(yGridPts.size(), 0.0));
	float buffer;

	for (size_t ix=0; ix<xGridPts.size(); ix++) {
		for (size_t iy=0; iy<yGridPts.size(); iy++) {
			file_in.read((char*)&buffer, sizeof(buffer));
            if (file_in.eof()) break;
			edens[ix][iy] = EoS.interp(static_cast<double>(buffer));
		}
	}

	file_in.close();

	edensint.SetData(xGridPts, yGridPts, edens);

	return 1;
}

static double CalcNumerator(interpFun &edensint, int n, double psin)
{
    double numerator = 0.0L;

	double x, y, rho, phi;

	for (size_t iq=0; iq<QMCPtsN; iq++)
	{
		x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[0][iq];
		y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[0][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		numerator += std::pow(rho, static_cast<double>(m))*std::cos(static_cast<double>(n)*(phi - psin))*edensint.interp(x, y);
	}

	numerator *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));

    return numerator;
}

static double CalcDenominator(interpFun &edensint)
{
    double denominator = 0.0L;

	double x, y, rho;

	for (size_t iq=0; iq<QMCPtsN; iq++)
	{
		x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[1][iq];
		y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[1][iq];
		rho = std::sqrt(x*x + y*y);
		denominator += std::pow(rho, static_cast<double>(m))*edensint.interp(x, y);
	}

	denominator *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));

    return denominator;
}

static std::map<size_t, std::map<int, double>> Epsn;

int CalcEpsn()
{
    #pragma omp parallel for
    for (size_t eventID=1; eventID<=eventN; eventID++)
    {
        interpFun edensInt; LoadEvol(eventID, edensInt);
		double denominator = CalcDenominator(edensInt);
        for (const auto& n : nList) {
			double numerator = CalcNumerator(edensInt, n, Psin[eventID][n]);
            Epsn[eventID][n] = -1.0*numerator/denominator;
        }
    }

    return 1;
}

int ExportEpsn()
{
	std::string path_out = "epsn_m=" + std::to_string(m) + ".dat";
    std::ofstream file_out(path_out, std::ios_base::out);
    if (!file_out.is_open()) {
        std::cerr << "Error: unable to open Epsn output file. Aborting..." << std::endl;
        return -1;
    }

    file_out << std::fixed << "#eventID ";
    for (const auto& n : nList)
        file_out << std::setw(11) << "Eps_" << n << " ";
    file_out << std::endl;

    for (size_t eventID=1; eventID<=eventN; eventID++) {
        file_out << std::fixed << std::setw(8) << eventID << " ";
        for (const auto& n : nList) {
            file_out << std::scientific << std::setw(11) << std::setprecision(6) << Epsn[eventID][n] << " ";
        }
        file_out << "\n";
    }

    file_out.close();

    return 1;
}