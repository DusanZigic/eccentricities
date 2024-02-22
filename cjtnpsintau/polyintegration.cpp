#include "polyintegration.hpp"

#include <vector>
#include <cmath>

static size_t LocatePoint(const std::vector<double> &data, double x, int interpolationOrder)
{
	size_t ju, jm, jl;
	size_t mm = interpolationOrder + 1;
	size_t n = data.size();
	size_t mmmin;
	bool ascnd = (data.back() >= data.front());
	jl = 0;
	ju = n - 1;
	while (ju - jl > 1)
	{
		jm = (ju + jl) >> 1;
		if ((x >= data[jm]) == ascnd) {
			jl = jm;
		}
		else {
			ju = jm;
		}
	}
	mmmin = n - mm < jl - ((mm - 2) >> 1) ? n - mm : jl - ((mm - 2) >> 1);
	return (0 > mmmin ? 0 : mmmin);
}

static void PolynomialCoeff(const std::vector<double> &dataX, const std::vector<double> &dataF, std::vector<double> &coeff)
{
	size_t n = dataX.size();

	coeff.resize(n); fill(coeff.begin(), coeff.end(), 0.0);

	size_t k, j, i;
	double phi, ff, b;
	
	std::vector<double> s(n); fill(s.begin(), s.end(), 0.0);
	
	s[n - 1] = -dataX[0];

	for (i = 1; i<n; i++)
	{
		for (j = n - 1 - i; j<n - 1; j++)
		{
			s[j] -= dataX[i] * s[j + 1];
		}
		s[n - 1] -= dataX[i];
	}

	for (j = 0; j<n; j++)
	{
		phi = n;
		
		for (k = n - 1; k>0; k--)
		{
			phi = k * s[k] + dataX[j] * phi;
		}
		
		ff = dataF[j] / phi;
		
		b = 1.0;
		
		for (k = n - 1; k >= 0; k--)
		{
			coeff[k] += b * ff;
			b = s[k] + dataX[j] * b;
		}
	}
}

double LinearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata)
{
	if (xdata.size() < 2) return 0.0;

	std::vector<double> k, c;
	for (size_t i=0; i<(xdata.size()-1); i++)
	{
		k.push_back((fdata[i+1]-fdata[i])/(xdata[i+1]-xdata[i]));
		c.push_back(fdata[i]-k.back()*xdata[i]);
	}

	double res = 0.0;

	for (size_t i=0; i<(xdata.size()-1); i++)
		res += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	return res;
}

double LinearIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit)
{
	if (xdata.size() < 2) return 0.0;

	std::vector<double> k, c;
	for (size_t i=0; i<(xdata.size()-1); i++)
	{
		k.push_back((fdata[i+1]-fdata[i])/(xdata[i+1]-xdata[i]));
		c.push_back(fdata[i]-k.back()*xdata[i]);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0L;

	for (size_t i=0; i<(xdata.size()-1); i++)
		sum += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	size_t lowLimitPos = LocatePoint(xdata, lowLimit, 1);

	double lowSum = 0.0L;

	for (size_t i=0; i<lowLimitPos; i++)
		lowSum += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);

	lowSum += 0.5*k[lowLimitPos]*(lowLimit*lowLimit - xdata[lowLimitPos]*xdata[lowLimitPos]) + c[lowLimitPos]*(lowLimit - xdata[lowLimitPos]);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	size_t highLimitPos = LocatePoint(xdata, highLimit, 1);

	double highSum = 0.0L;

	highSum += 0.5*k[highLimitPos]*(xdata[highLimitPos+1]*xdata[highLimitPos+1] - highLimit*highLimit) + c[highLimitPos]*(xdata[highLimitPos+1] - highLimit);

	for (size_t i=highLimitPos+1; i<xdata.size()-1; i++)
		highSum += 0.5*k[i]*(xdata[i+1]*xdata[i+1] - xdata[i]*xdata[i]) + c[i]*(xdata[i+1] - xdata[i]);	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - highSum - lowSum);
}

double CubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata)
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	std::vector<std::vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (size_t i=0; i<coefficents.size(); i++)
	{
		size_t pointLocation = LocatePoint(xdata, xdata[i], 3);
		std::vector<double> xdatatemp(xdata.begin()+pointLocation, xdata.begin()+pointLocation+4);
		std::vector<double> fdatatemp(fdata.begin()+pointLocation, fdata.begin()+pointLocation+4);
		PolynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//calculating value of integral:
	double sum = 0.0L;
	for (size_t i=0; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			sum += 1.0/(j+1)*coefficents[i][j]*(std::pow(xdata[i+1], static_cast<double>(j+1)) - std::pow(xdata[i], static_cast<double>(j+1)));

	return sum;
}

double CubicIntegrate(const std::vector<double> &xdata, const std::vector<double> &fdata, double lowLimit, double highLimit)
{
	if (xdata.size() < 2) return 0;

	//calculating polynomial coefficients for each segment:
	std::vector<std::vector<double>> coefficents; coefficents.resize(xdata.size() - 1);
	for (size_t i=0; i<coefficents.size(); i++)
	{
		size_t pointLocation = LocatePoint(xdata, xdata[i], 3);
		std::vector<double> xdatatemp(xdata.begin()+pointLocation, xdata.begin()+pointLocation+4);
		std::vector<double> fdatatemp(fdata.begin()+pointLocation, fdata.begin()+pointLocation+4);
		PolynomialCoeff(xdatatemp, fdatatemp, coefficents[i]);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of full integral (in it's whole range):
	double sum = 0.0L;
	for (size_t i=0; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			sum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                        (std::pow(xdata[i+1], static_cast<double>(j+1)) - std::pow(xdata[i], static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from lower range to lower limit:
	size_t lowLimitPos = LocatePoint(xdata, lowLimit, 1);

	double lowSum = 0.0L;

	for (size_t i=0; i<lowLimitPos; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			lowSum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                            (std::pow(xdata[i+1], static_cast<double>(j+1)) - std::pow(xdata[i], static_cast<double>(j+1)));

	for (size_t j=0; j<coefficents[lowLimitPos].size(); j++)
		lowSum += 1.0L/static_cast<double>(j+1)*coefficents[lowLimitPos][j]*
                        (std::pow(lowLimit, static_cast<double>(j+1)) - std::pow(xdata[lowLimitPos], static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//calculating value of integral from higher limit to higer range:
	size_t highLimitPos = LocatePoint(xdata, highLimit, 1);

	double highSum = 0.0L;

	for (size_t j=0; j<coefficents[highLimitPos].size(); j++)
		highSum += 1.0L/static_cast<double>(j+1)*coefficents[highLimitPos][j]*
                        (std::pow(xdata[highLimitPos+1], static_cast<double>(j+1)) - std::pow(highLimit, static_cast<double>(j+1)));

	for (size_t i=highLimitPos+1; i<xdata.size()-1; i++)
		for (size_t j=0; j<coefficents[i].size(); j++)
			highSum += 1.0L/static_cast<double>(j+1)*coefficents[i][j]*
                            (std::pow(xdata[i+1], static_cast<double>(j+1)) - std::pow(xdata[i], static_cast<double>(j+1)));

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//integral value is full-low-high
	return (sum - highSum - lowSum);
}