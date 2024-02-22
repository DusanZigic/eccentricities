#include "arsenal.hpp"
#include "polyintegration.hpp"
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

static const double TCRIT = 0.155L;

static size_t eventN; static std::vector<unsigned int> nList; static unsigned int m; size_t QMCPtsN;

int GetInputs(int argc, char const *argv[])
{
    std::vector<std::string> inputs; for (int i=1; i<argc; i++) inputs.push_back(argv[i]);

	if ((inputs.size() == 1) && ((inputs[0] == "-h") || (inputs[0] == "--h"))) {
		std::cout << "default values: --eventN=1000 --n=1-8 --m=2 --QMCPtsN=5000000" << std::endl;
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
	std::string nstr =   "1-8"; if (inputparams.count("n")       > 0)    nstr = 		  inputparams["n"];
    			   m =       2; if (inputparams.count("m")       > 0)       m = std::stoi(inputparams["m"]);
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

static size_t phiGaussPtsN;
static std::map<std::string, std::vector<double>> phiGausPts;

int LoadPhiGaussPts()
{
    phiGausPts["phi"]     = std::vector<double>(0);
    phiGausPts["weights"] = std::vector<double>(0);

    std::ifstream file_in("./phigausspts.dat", std::ios_base::in);
    if(!file_in.is_open()) {
        std::cerr << "Error: unable to open phi Gauss points file. Aborting..." << std::endl;
        return -1;
    }

    std::string line; double buffer;

    while(std::getline(file_in, line))
    {
        if (line.at(0) == '#')
            continue;
        
        std::stringstream ss(line);
        ss >> buffer; phiGausPts["phi"].push_back(buffer);
        ss >> buffer; phiGausPts["weights"].push_back(buffer);
    }

    phiGaussPtsN = phiGausPts["phi"].size();

    file_in.close();

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

    for (int i_tau=0; i_tau<tau_step_n; i_tau++) tauGridPts.push_back(tau_min + static_cast<double>(i_tau)*tau_step);
    for (int i_x=0;   i_x<x_step_n;     i_x++)     xGridPts.push_back(  x_min + static_cast<double>(  i_x)*  x_step);
    for (int i_y=0;   i_y<y_step_n;     i_y++)     yGridPts.push_back(  y_min + static_cast<double>(  i_y)*  y_step);

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

static int LoadEvol(size_t event_id, std::vector<double> &taupts, interpFun &tempsint, interpFun &edensint)
{
    std::string path_in = "./tempevols/TProfile_" + std::to_string(event_id) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in | std::ios_base::binary);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open evolution file for event " + std::to_string(event_id) + ". Aborting..." << std::endl;
        return -1;
    }

    std::vector<double> temps, edens;
    float buffer;

    while (true)
    {
        file_in.read((char*)&buffer, sizeof(buffer));
        if (file_in.eof()) break;
        temps.push_back(static_cast<double>(buffer));
        edens.push_back(EoS.interp(static_cast<double>(buffer)));
    }

    file_in.close();

    tempsint.SetData(tauInterpPts, xInterpPts, yInterpPts, temps);
    edensint.SetData(tauInterpPts, xInterpPts, yInterpPts, edens);

    double maxTemp = temps[0], x0 = 0.0L, y0 = 0.0L;
    for (const auto &x : xGridPts) {
        for (const auto &y : yGridPts) {
            double temp = tempsint.interp(tauGridPts[0], x, y);
            if (temp > maxTemp) {
                maxTemp = temp;
                x0 = x;
                y0 = y;
            }
        }
    }

    size_t tauCnt = 0;
    while (tempsint.interp(tauGridPts[tauCnt], x0, y0) > TCRIT) {
        taupts.push_back(tauGridPts[tauCnt]);
        tauCnt++;
    }    

    return 1;
}

static int LoadBCPoints(size_t event_id, std::vector<std::vector<double>> &bcpts)
{
    std::string path_in = "./bcpoints/BinaryCollPoints_" + std::to_string(event_id) + ".dat";
    std::ifstream file_in(path_in, std::ios_base::in);
    if (!file_in.is_open()) {
        std::cerr << "Error: unable to open bc points file for event " + std::to_string(event_id) + ". Aborting..." << std::endl;
        return -1;
    }

    std::string line; double bufferX, bufferY;

    while(std::getline(file_in, line))
    {
        if (line.at(0) == '#')
            continue;
        
        std::stringstream ss(line);
        ss >> bufferX; ss >> bufferY;
        bcpts.push_back(std::vector<double>{bufferX, bufferY});
    }

    file_in.close();

    return 1;
}

static int CalcjTnTauPhi(const std::vector<double> &taupts, interpFun &tempsint, const std::vector<std::vector<double>> &bcpoints, std::vector<std::vector<double>> &jtntauphi)
{
    jtntauphi.resize(taupts.size(), std::vector<double>(phiGaussPtsN, 0.0L));

    double x, y, tau, phi; size_t Nbc; double temp;
    const std::vector<double> xDomain = tempsint.domain()[1], yDomain = tempsint.domain()[2];

    for (size_t itau=0; itau<taupts.size(); itau++) {
        tau = taupts[itau];
        for (size_t iphi=0; iphi<phiGaussPtsN; iphi++) {
            phi = phiGausPts["phi"][iphi];
            Nbc = 0;
            for (const auto& bcpoint : bcpoints) {
                x = bcpoint[0] + tau*std::cos(phi);
                y = bcpoint[1] + tau*std::sin(phi);
                if ((x >= xDomain[0]) && (x <= xDomain[1]) && (y >= yDomain[0]) && (y <= yDomain[1])) {
                    Nbc++;
                    temp = tempsint.interp(tau, x, y);
                    jtntauphi[itau][iphi] += temp*temp*temp;
                }
            }
            jtntauphi[itau][iphi] /= static_cast<double>(Nbc);
        }
    }

    return 1;
}

static int CalcPsin(const std::vector<double> &taupts, interpFun &edensint, unsigned int n, std::vector<double> &psintau)
{
    for (const auto &tau : taupts) {
        double sinsum = 0.0L, cossum = 0.0L;
        double x, y, rho, phi;
        for (size_t iq=0; iq<QMCPtsN; iq++) {
            x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[0][iq];
            y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[0][iq];
            rho = std::sqrt(x*x + y*y);
            phi = std::atan2(y, x);
            sinsum += std::pow(rho, static_cast<double>(m))*std::sin(static_cast<double>(n)*phi)*edensint.interp(tau, x, y);

            x = xGridPts.front() + (xGridPts.back() - xGridPts.front())*QMCPtsX[1][iq];
            y = yGridPts.front() + (yGridPts.back() - yGridPts.front())*QMCPtsY[1][iq];
            rho = std::sqrt(x*x + y*y);
            phi = std::atan2(y, x);
            cossum += std::pow(rho, static_cast<double>(m))*std::cos(static_cast<double>(n)*phi)*edensint.interp(tau, x, y);
        }
        sinsum *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));
        cossum *= ((xGridPts.back() - xGridPts.front())*(yGridPts.back() - yGridPts.front())/static_cast<double>(QMCPtsN));
        psintau.push_back(1.0L/static_cast<double>(n)*std::atan2(sinsum, cossum) + PI/static_cast<double>(n));
    }

    return 1;
}

static int CalcjTnTau(const std::vector<double> &taupts, std::vector<std::vector<double>> &jtntauphi, unsigned int n, std::vector<double> &psintau, std::vector<double> &jtntau)
{
    jtntau = std::vector<double>(taupts.size(), 0.0L);
    for (size_t itau=0; itau<taupts.size(); itau++) {
        double numerator = 0.0L, denominator = 0.0L;
        for (size_t iphi=0; iphi<phiGaussPtsN; iphi++) {
            numerator   += phiGausPts["weights"][iphi]*std::cos(static_cast<double>(n)*(phiGausPts["phi"][iphi] - psintau[itau]))*jtntauphi[itau][iphi];
            denominator += phiGausPts["weights"][iphi]*jtntauphi[itau][iphi];
        }
        jtntau[itau] = -1.0L*numerator/denominator;
    }

    return 1;
}

static std::map<size_t, std::map<unsigned int, double>> jTn;

int CalcjTn()
{
    #pragma omp parallel for
    for (size_t eventID=1; eventID<=eventN; eventID++)
    {
        std::vector<double> tauPts; interpFun tempsInt, edensInt; LoadEvol(eventID, tauPts, tempsInt, edensInt);
        std::vector<std::vector<double>> bcPoints;                LoadBCPoints(eventID, bcPoints);
        std::vector<std::vector<double>> jTnTauPhi;               CalcjTnTauPhi(tauPts, tempsInt, bcPoints, jTnTauPhi);
        for (const auto &n : nList) {
            std::vector<double> PsinTau; CalcPsin(tauPts, edensInt, n, PsinTau);
            std::vector<double> jTnTau;  CalcjTnTau(tauPts, jTnTauPhi, n, PsinTau, jTnTau);
            jTn[eventID][n] = 1.0L/(tauPts.back() - tauPts.front())*LinearIntegrate(tauPts, jTnTau);
        }
    }

    return 1;
}

int ExportjTn()
{
    std::ofstream file_out("./jTnQMC.dat", std::ios_base::out);
    if (!file_out.is_open()) {
        std::cerr << "Error: unable to open jTn output file. Aborting..." << std::endl;
        return -1;
    }

    file_out << std::fixed << "#eventID ";
    for (const auto& n : nList)
        file_out << std::setw(12) << "jT_" << n << " ";
    file_out << std::endl;

    for (size_t eventID=1; eventID<=eventN; eventID++) {
        file_out << std::fixed << std::setw(8) << eventID << " ";
        for (const auto& n : nList) {
            file_out << std::scientific << std::setw(13) << std::setprecision(6) << jTn[eventID][n] << " ";
        }
        file_out << "\n";
    }

    file_out.close();
    
    return 1;
}