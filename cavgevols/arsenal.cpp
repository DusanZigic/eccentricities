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

static std::string rotationAngle;
static std::vector<std::string> rotationAngleList;
static size_t eventN;

int GetInputs(int argc, char const *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && ((inputs[0] == "-h") || (inputs[0] == "--h"))) {
		std::cout << "default values: --rotationAngle=plus-angle-minus-pi-n-pi-2" << std::endl;
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

    rotationAngle = "plus-angle-minus-pi-n-minus-pi-2"; if (inputparams.count("rotationAngle") > 0) rotationAngle = inputparams["rotationAngle"];

    rotationAngleList = {"minus-angle", "minus-angle-minus-pi-n", "minus-angle-plus-pi-n",
                          "plus-angle",  "plus-angle-minus-pi-n",  "plus-angle-plus-pi-n"};
    size_t rotationAngleListSize = rotationAngleList.size();
    for (size_t i=0; i<rotationAngleListSize; i++) rotationAngleList.push_back(rotationAngleList[i] + "-minus-pi-2");

    if (std::find(rotationAngleList.begin(), rotationAngleList.end(), rotationAngle) == rotationAngleList.end()) {
        std::cerr << "Error: rotation angle not an available option. Aborting..." << std::endl;
        return -1;
    }

    eventN = 0; if (inputparams.count("eventN") > 0) eventN = std::stoi(inputparams["eventN"]);

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

static std::map<size_t, std::map<int, double>> Psin;
static std::vector<int> nList;

int LoadPsin()
{
    std::string path_in = "./psin.dat";
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

    if (eventN == 0)
        eventN = Psin.size();

    return 1;
}

static size_t minEvolLength;

int GetMinEvolLength()
{
    std::vector<size_t> evolLengths;

    for (size_t eventID=1; eventID<=eventN; eventID++)
    {
        std::string path_in = "./tempevols/TProfile_" + std::to_string(eventID) + ".dat";
        std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
        if (!file_in.is_open()) {
            std::cerr << "Error: unable to open evolution file for event " + std::to_string(eventID) + ". Aborting..." << std::endl;
            return -1;
        }
        file_in.seekg(0, std::ios_base::end);
        evolLengths.push_back(file_in.tellg() * 8 / 32); //size in bytes * 8bits / 32bit float
        file_in.close();
    }

    minEvolLength = *std::min_element(evolLengths.begin(), evolLengths.end());

    return 1; 
}

static std::vector<double> tauInterpPts, xInterpPts, yInterpPts;
static std::vector<double> tauGridPts,   xGridPts,   yGridPts;

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
		double tau = tau_min; while (tau <= tau_max) {tau_step_n++; tau+=tau_step;}
	}

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buff; x_min  = buff;
							   		ss >> buff; x_max  = buff;
							   		ss >> buff; x_step = buff;
		double x = x_min; while (x <= x_max) {x_step_n++; x+=x_step;}
	}

	{
		std::getline(file_in, line);
		std::getline(file_in, line);
		std::stringstream ss(line); ss >> buff; y_min  = buff;
							   		ss >> buff; y_max  = buff;
							   		ss >> buff; y_step = buff;
		double y = y_min; while (y <= y_max) {y_step_n++; y+=y_step;}
	}

	file_in.close();

    for (int i_tau=0; i_tau<tau_step_n; i_tau++) {
	    for (int i_x=0; i_x<x_step_n; i_x++) {
	        for (int i_y=0; i_y<y_step_n; i_y++) {
                tauInterpPts.push_back(tau_min + static_cast<double>(i_tau)*tau_step);
                  xInterpPts.push_back(  x_min + static_cast<double>(  i_x)*  x_step);
                  yInterpPts.push_back(  y_min + static_cast<double>(  i_y)*  y_step);
            }
        }
    }

    tauInterpPts.resize(minEvolLength);
      xInterpPts.resize(minEvolLength);
      yInterpPts.resize(minEvolLength);

    for (int i_tau=0; i_tau<tau_step_n; i_tau++) tauGridPts.push_back(tau_min + static_cast<double>(i_tau)*tau_step);
    for (int i_x=0;   i_x<x_step_n;     i_x++)     xGridPts.push_back(  x_min + static_cast<double>(  i_x)*  x_step);
    for (int i_y=0;   i_y<y_step_n;     i_y++)     yGridPts.push_back(  y_min + static_cast<double>(  i_y)*  y_step);

    tauGridPts.resize(minEvolLength / xGridPts.size() / yGridPts.size());

	return 1;
}

static int LoadEvol(size_t event_id, interpFun& tempsint, interpFun& edensint)
{
    std::string path_in = "./tempevols/TProfile_" + std::to_string(event_id) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open evolution file for event " + std::to_string(event_id) + ". Aborting..." << std::endl;
        return -1;
    }

    std::vector<double> temps, edens;
    float buffer;
    size_t counter = 0;

    while (true)
    {
        file_in.read((char*)&buffer, sizeof(buffer));
        if (file_in.eof()) break;
        temps.push_back(static_cast<double>(buffer));
        edens.push_back(EoS.interp(static_cast<double>(buffer)));
        counter++;
        if (counter >= minEvolLength) break;
    }

    file_in.close();

    tempsint.SetData(tauInterpPts, xInterpPts, yInterpPts, temps);
    edensint.SetData(tauInterpPts, xInterpPts, yInterpPts, edens);

    return 1;
}

static int ExportRotatedEvols(int n, std::vector<std::vector<std::vector<double>>> &avgtemps, std::vector<std::vector<std::vector<double>>> &avgedens)
{
    std::ofstream file_out("avgrotevoln" + std::to_string(n) + ".dat");
    if (!file_out.is_open()) {
        std::cerr << "Error: unable to open evolutions output file. Aborting..." << std::endl;
        return -1;
    }

    file_out << "#";
    file_out << std::fixed << setw(5)  <<   "tau" << " ";
    file_out << std::fixed << setw(5)  <<     "x" << " ";
    file_out << std::fixed << setw(5)  <<     "y" << " ";
    file_out << std::fixed << setw(12) <<  "temp" << " ";
    file_out << std::fixed << setw(12) << "edens" << "\n";

    for (size_t itau=0; itau<tauGridPts.size(); itau++) {
        for (size_t ix=0; ix<xGridPts.size(); ix++) {
            for (size_t iy=0; iy<xGridPts.size(); iy++) {
                file_out << std::fixed      << std::setw(6)  << std::setprecision(3) <<       tauGridPts[itau] << " ";
                file_out << std::fixed      << std::setw(5)  << std::setprecision(1) <<           xGridPts[ix] << " ";
                file_out << std::fixed      << std::setw(5)  << std::setprecision(1) <<           yGridPts[iy] << " ";
                file_out << std::scientific << std::setw(12) << std::setprecision(6) << avgtemps[itau][ix][iy] << " ";
                file_out << std::scientific << std::setw(12) << std::setprecision(6) << avgedens[itau][ix][iy] << "\n";
            }
        }
    }

    file_out.close();

    return 1;
}

int AvgRotatedEvols()
{
    for (const auto& n : nList)
    {
        std::vector<std::vector<std::vector<double>>> avgTemps(tauGridPts.size(), std::vector<vector<double>>(xGridPts.size(), std::vector<double>(yGridPts.size(), 0.0L)));
        std::vector<std::vector<std::vector<double>>> avgEdens(tauGridPts.size(), std::vector<vector<double>>(xGridPts.size(), std::vector<double>(yGridPts.size(), 0.0L)));

        for (size_t eventID=1; eventID<=eventN; eventID++)
        {
            interpFun tempsInt, edensInt; if (LoadEvol(eventID, tempsInt, edensInt) != 1) return -1;

            double tau, x, y, rotAngle = 0.0;

            if (rotationAngle == rotationAngleList[0])  rotAngle = -1.0L*Psin[eventID][n];
            if (rotationAngle == rotationAngleList[1])  rotAngle = -1.0L*Psin[eventID][n] - PI/static_cast<double>(n);
            if (rotationAngle == rotationAngleList[2])  rotAngle = -1.0L*Psin[eventID][n] + PI/static_cast<double>(n);
            if (rotationAngle == rotationAngleList[3])  rotAngle =       Psin[eventID][n];
            if (rotationAngle == rotationAngleList[4])  rotAngle =       Psin[eventID][n] - PI/static_cast<double>(n);
            if (rotationAngle == rotationAngleList[5])  rotAngle =       Psin[eventID][n] + PI/static_cast<double>(n);
            if (rotationAngle == rotationAngleList[6])  rotAngle = -1.0L*Psin[eventID][n]                            - PI/2.0L;
            if (rotationAngle == rotationAngleList[7])  rotAngle = -1.0L*Psin[eventID][n] - PI/static_cast<double>(n)- PI/2.0L;
            if (rotationAngle == rotationAngleList[8])  rotAngle = -1.0L*Psin[eventID][n] + PI/static_cast<double>(n)- PI/2.0L;
            if (rotationAngle == rotationAngleList[9])  rotAngle =       Psin[eventID][n]                            - PI/2.0L;
            if (rotationAngle == rotationAngleList[10]) rotAngle =       Psin[eventID][n] - PI/static_cast<double>(n)- PI/2.0L;
            if (rotationAngle == rotationAngleList[11]) rotAngle =       Psin[eventID][n] + PI/static_cast<double>(n)- PI/2.0L;

            for (size_t itau=0; itau<tauGridPts.size(); itau++) {
                tau = tauGridPts[itau];
                for (size_t ix=0; ix<xGridPts.size(); ix++) {
                    for (size_t iy=0; iy<xGridPts.size(); iy++) {
                        x = xGridPts[ix]*std::cos(rotAngle) - yGridPts[iy]*std::sin(rotAngle);
                        y = xGridPts[ix]*std::sin(rotAngle) + yGridPts[iy]*std::cos(rotAngle);
                        if ((x >= xGridPts.front()) && (x <= xGridPts.back()) && (y >= yGridPts.front()) && (y <= yGridPts.back())) {
                            avgTemps[itau][ix][iy] += tempsInt.interp(tau, x, y);
                            avgEdens[itau][ix][iy] += edensInt.interp(tau, x, y);
                        }
                    }
                }
            }
        }

        for (size_t itau=0; itau<tauGridPts.size(); itau++) {
            for (size_t ix=0; ix<xGridPts.size(); ix++) {
                for (size_t iy=0; iy<xGridPts.size(); iy++) {
                    avgTemps[itau][ix][iy] /= static_cast<double>(eventN);
                    avgEdens[itau][ix][iy] /= static_cast<double>(eventN);
                }
            }
        }

        if (ExportRotatedEvols(n, avgTemps, avgEdens) != 1) return -2;
    }

    return 1;
}