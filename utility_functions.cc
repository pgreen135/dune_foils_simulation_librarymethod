#include "utility_functions.h"
#include <cmath>
#include "TMath.h"
#include "TVector3.h"
#include "TF1.h"
#include "TRandom.h"
#include <chrono>
#include <iostream>

//Poisson Distribution
int utility::poisson(double mean, double draw, double eng)
{
	int number = 0;
	const int border = 16;
	double limit = 2e9;

	if(mean <= border) {
		double position = draw;
		double poissonValue = std::exp(-mean);
		double poissonSum = poissonValue;

		while(poissonSum <= position) {
			number++;
			poissonValue *= mean/number;
			poissonSum += poissonValue;
		}
		return number;
	} // the case of mean <= 16

	double value, t, y;
	t = std::sqrt(-2*std::log(draw));
	y = 2*3.141592654*eng;
	t *= std::cos(y);
	value = mean + t*std::sqrt(mean) + 0.5;
	if(value <= 0) {return 0; }
	if(value >= limit) { return limit; }
	return value;
}


//Beta decay function
double utility::SpectrumFunction(double *x, double *par)
{
	double KE = *x;
	double Q = par[0];
	double MassE = 0.510998910; // mass electron - MeV/c^2

	double N = std::sqrt(pow(KE,2) + 2* KE * MassE ) * std::pow( (Q - KE), 2 ) * (KE + MassE );

	return N;
}

//Supernova neutrino energy spectrum
//Spectrum energy comes from a 1 kpc galactic supernova
double utility::fsn(double *x, double *par)
{

	double E = *x;
	double Eav = par[0];

	double f_nu = std::pow(E,3)*std::exp(-4*E/Eav);

	return f_nu;

}

//Solar neutrino energy spectrum
double utility::fso(double *x, double *par)
{

        double E = *x;
        double Eav = par[0];

        double f_nu = -0.00146609+0.00589747*E+0.0108878*std::pow(E,2)-0.00174571*std::pow(E,3)+4.97905e-05*std::pow(E,4)+3.15376e-06*std::pow(E,5)-1.37052e-07*std::pow(E,6);//pol6 fit with chi square = 0.00013515

        return std::abs(f_nu);

}

//Radon-222 decay energy spectrum (Gaussian about the alpha decay energy)

double utility::Rn_function(double *x, double *par)
{

	double E = *x;
	double Q_Rn = par[0];

	double sigma = 0.01;
	//double sigma = (E - Q_Rn) / std::sqrt(1.3863); // 1,3863 = ln(4)
	double sigma_sq = sigma * sigma;

	double gauss = 1/(sigma * std::sqrt(2*3.1416)) * std::exp((-1*std::pow((E-Q_Rn),2))/(2*sigma_sq));

	return gauss;

}


//Scintillation decay function
double utility::Scintillation_function(double *t, double *par){

	double time = *t;
	double t_singlet = par[0];
	double t_triplet = par[1];
	double type = par[2]; // type will be defined at 0 or 1, 0 is an electron, 1 is an alpha particle
	double singlet_part;
	double triplet_part;

	if(type == 0){ // particle is an electron
	  singlet_part = 0.25;
	  triplet_part = 0.75;
	}

	if(type == 1){ // particle is an alpha
	  singlet_part = 0.75;
	  triplet_part = 0.25;
	}

	double Scintillation = exp(-(time/t_singlet))*singlet_part/t_singlet;  + exp(-(time/t_triplet))*triplet_part/t_triplet;

	return Scintillation;

}