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

static size_t eventN; static unsigned int m; static std::vector<int> nList; static size_t QMCPtsN;

int GetInputs(int argc, char const *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && ((inputs[0] == "-h") || (inputs[0] == "--h"))) {
		std::cout << "default values: --eventN=1000 --m=2 --n=1-8 --QMCPtsN=5000000" << std::endl;
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

    		  eventN =    1000; if (inputparams.count("eventN")  > 0)  eventN = std::stoi(inputparams["eventN"]);
    			   m =       2; if (inputparams.count("m")       > 0)       m = std::stoi(inputparams["m"]);
	std::string nstr =   "1-8"; if (inputparams.count("n")       > 0)    nstr = 		  inputparams["n"];
    		 QMCPtsN = 5000000; if (inputparams.count("QMCPtsN") > 0) QMCPtsN = std::stoi(inputparams["QMCPtsN"]);

    if (eventN < 0) {
		std::cerr << "Error: provided eventN parameter must be positive integer. Aborting...\n";
		return -1;
	}

	std::vector<int> nRange;
	if (nstr.find("-") != std::string::npos) {
		std::stringstream nsstr(nstr);
		for (int n; nsstr >> n;) {
        	nRange.push_back(n);
        	if (nsstr.peek() == '-')
            	nsstr.ignore();
    	}
	}
	else {
		std::cerr << "Error: numbers in n parameter must be separated by: '-'. Aborting..." << std::endl;
		return -2;
	}

	for (int n=nRange.front(); n<=nRange.back(); n++)
		nList.push_back(n);

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

static int CalcPsin(interpFun& edensint, int n, double& psin)
{
	double sinsum = 0.0L, cossum = 0.0L;

	double x, y, rho, phi;

	for (size_t iq=0; iq<QMCPtsN; iq++)
	{
		x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[0][iq];
		y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[0][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		sinsum += std::pow(rho, static_cast<double>(m))*std::sin(static_cast<double>(n)*phi)*edensint.interp(x, y);

		x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[1][iq];
		y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[1][iq];
		rho = std::sqrt(x*x + y*y);
		phi = std::atan2(y, x);
		cossum += std::pow(rho, static_cast<double>(m))*std::cos(static_cast<double>(n)*phi)*edensint.interp(x, y);
	}

	sinsum *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));
	cossum *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));

	psin = 1.0L/static_cast<double>(n) * std::atan2(sinsum, cossum) + PI/static_cast<double>(n);

	return 1;
}

static std::map<size_t, std::map<int, double>> Psin;

int CalcPsin()
{
    #pragma omp parallel for
    for (size_t eventID=1; eventID<=eventN; eventID++)
    {
        interpFun edensInt; LoadEvol(eventID, edensInt);
        for (const auto& n : nList) {
            CalcPsin(edensInt, n, Psin[eventID][n]);
        }
    }

	return 1;
}

int ExportPsin()
{
    std::ofstream file_out("psin.dat", std::ios_base::out);
    if (!file_out.is_open()) {
        std::cerr << "Error: unable to open Psin output file. Aborting..." << std::endl;
        return -1;
    }

    file_out << std::fixed << "#eventID ";
    for (const auto& n : nList)
        file_out << std::setw(11) << "Psi_" << n << " ";
    file_out << std::endl;

    for (size_t eventID=1; eventID<=eventN; eventID++) {
        file_out << std::fixed << std::setw(8) << eventID << " ";
        for (const auto& n : nList) {
            file_out << std::scientific << std::setw(11) << std::setprecision(6) << Psin[eventID][n] << " ";
        }
        file_out << "\n";
    }

    file_out.close();

    return 1;
}