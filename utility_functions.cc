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


/*

------------------ not currently in use ---------------------


//        ======TIMING PARAMETRIZATION=====
//For more info, Diego wrote these functions.
double utility::finter_d(double *x, double *par) {

	double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
	double y2 = TMath::Exp(par[3]+x[0]*par[4]);

	return TMath::Abs(y1 - y2);
}

double utility::LandauPlusExpoFinal(double *x, double *par)
{
	// par0 = joining point
	// par1 = Landau MPV
	// par2 = Landau width
	// par3 = normalization
	// par4 = Expo cte
	// par5 = Expo tau

	double y1 = TMath::Landau(x[0],par[1],par[2])*par[3];
	double y2 = TMath::Exp(par[4]+x[0]*par[5]);
	if(x[0] > par[0]) y1 = 0.;
	if(x[0] < par[0]) y2 = 0.;

	return (y1 + y2);
}

double utility::finter_r(double *x, double *par) {

	double y1 = par[2]*TMath::Landau(x[0],par[0],par[1]);
	double y2 = par[5]*TMath::Landau(x[0],par[3],par[4]);

	return TMath::Abs(y1 - y2);
}

double utility::LandauPlusLandauFinal(double *x, double *par)
{
	// par0 = joining point
	// par1 = Landau1 MPV
	// par2 = Landau1 width
	// par3 = Landau1 normalization
	// par4 = Landau2 MPV
	// par5 = Landau2 width
	// par6 = Landau2 normalization

	double y1 = par[3]*TMath::Landau(x[0],par[1],par[2]);
	double y2 = par[6]*TMath::Landau(x[0],par[4],par[5]);
	if(x[0] > par[0]) y1 = 0.;
	if(x[0] < par[0]) y2 = 0.;

	return (y1 + y2);
}
//-----------------------------------------------------------
// Parametrization of the VUV light timing (result from direct transport + Rayleigh scattering ONLY)
// using a landau + expo function.The function below returns the arrival time distribution given the
// distance IN CENTIMETERS between the scintillation/ionization point and the optical detector.
/*std::vector<double> utility::GetVUVTime(double distance, int number_photons) {
   std::vector<double> arrival_time_distrb;
   return arrival_time_distrb;
   }*/
/*
std::vector<double> utility::GetVUVTime(double distance, int number_photons) {
	//-----Distances in cm and times in ns-----//

	//std::cout << "Distance = " << distance << std::endl;

	//gRandom->SetSeed(0);
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
	std::vector<double> arrival_time_distrb;
	arrival_time_distrb.clear();
	arrival_time_distrb.reserve(number_photons);
	
	//range of distances where the parametrization is valid [~10 - 500cm], then:
    const double d_break = 500.;
    const double d_max = 750.;

	if(distance < 10 ) {
		arrival_time_distrb.assign(number_photons,10.);
		std::cout<<"WARNING: Parametrization failure - distance < 10 cm"<<std::endl;
		//std::cout<<"Too close/far to the PMT  -> set 0 VUV photons(?)!!!!!!"<<std::endl;
		return arrival_time_distrb;
	}
	if(distance > d_max) {
	    // TEMP FOR DUNE: If outside parametrization range, just set times at 100.		THIS NEEDS CHANGING 
	    arrival_time_distrb.assign(number_photons,10000.);
	    std::cout << "WARNING: Parametrization failure - distance > d_max" << std::endl;
	    return arrival_time_distrb;
	}
	
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
	// Defining the function for the different parameters:
    TF1 fparslogNorm ("fparslogNorm","pol7",0,d_break);
    TF1 fparsMPV ("fparsMPV","pol4",0,d_break);
    TF1 fparsWidth ("fparsWidth","pol3",0,d_break);
    TF1 fparsCte ("fparsCte","pol6",0,d_break);
    TF1 fparsSlope ("fparsSlope","pol7",0,d_break);

    // At long distances we extrapolate the behaviour of the parameters:
    TF1 fparslogNorm_far ("fparslogNorm_far","expo",d_break, d_max);
    double landauNormpars_far[2] = {2.23151, -0.00627503};
    TF1 fparsMPV_far ("fparsMPV_far","pol1",d_break, d_max);
    double landauMPVpars_far[2] = {-3.04952, 0.128638};
    TF1 fparsCte_far ("fparsCte_far","expo",d_break-50., d_max);
    double expoCtepars_far[2] = {3.69578, -0.00989582};
    // Parametrization data:
    double landauNormpars[8] = {7.85903, -0.108075, 0.00110999, -6.90009e-06,
		                        2.52576e-08, -5.39078e-11, 6.20863e-14, -2.97559e-17};
    double landauMPVpars[5] = {1.20259, 0.0582674, 0.000308053, -2.71782e-07, -3.37159e-10};
    double landauWidthpars[4] = {0.346667, -0.00768231, 0.000211825, -3.81361e-07};
    double expoCtepars[7] = {13.6592, -0.188798, 0.00192431, -1.10689e-05, 3.38425e-08,
		                     -5.20737e-11, 3.17657e-14};
    double expoSlopepars[8] = {-0.57011, 0.0156393, -0.000197461, 1.34491e-06, -5.24544e-09,
		                       1.1703e-11, -1.38811e-14, 6.78368e-18};

    // This ensures the parameter functions are setup only once ever
    bool setup_VUV_params = true;

    std::vector<std::vector<double>> VUVTRange = { { 0.30000000, 12.300000 }, { 1.0000000, 18.100000 }, { 1.7000000, 26.000000 }, { 2.2000000, 35.900000 }, { 3.0000000, 48.300000 }, { 3.1000000, 62.800000 }, { 4.0000000, 79.300000 }, { 4.0000000, 96.400000 }, { 4.1000000, 112.60000 }, { 4.3000000, 127.30000 }, { 5.0000000, 139.90000 }, { 5.1000000, 150.10000 }, { 5.4000000, 159.40000 }, { 5.9000000, 167.90000 }, { 5.2000000, 177.00000 }, { 6.0000000, 186.30000 }, { 5.0000000, 197.90000 }, { 6.4000000, 209.70000 }, { 6.5000000, 224.20000 }, { 7.1000000, 239.60000 }, { 6.8000000, 255.10000 }, { 7.1000000, 272.00000 }, { 7.5000000, 286.90000 }, { 7.7000000, 301.40000 }, { 8.1000000, 312.10000 }, { 8.9000000, 320.60000 }, { 10.000000, 325.60000 }, { 8.9000000, 328.00000 }, { 10.700000, 327.80000 }, { 11.400000, 326.90000 }, { 11.300000, 325.00000 }, { 11.300000, 324.80000 }, { 13.900000, 324.90000 }, { 14.800000, 327.40000 }, { 14.300000, 330.10000 }, { 16.200000, 336.80000 }, { 18.600000, 343.00000 }, { 19.000000, 350.20000 }, { 21.200000, 358.40000 }, { 22.800000, 362.40000 }, { 25.000000, 366.30000 }, { 26.900000, 366.50000 }, { 27.300000, 361.40000 }, { 30.400000, 352.90000 }, { 32.500000, 341.60000 }, { 35.100000, 345.40000 }, { 36.400000, 347.40000 }, { 40.300000, 348.60000 }, { 42.900000, 350.20000 }, { 45.400000, 351.90000 }, { 56.000000, 361.40000 }, { 57.100000, 362.90000 }, { 58.300000, 364.00000 }, { 60.000000, 365.50000 }, { 61.000000, 366.90000 }, { 62.100000, 368.20000 }, { 63.200000, 369.10000 }, { 64.500000, 370.60000 }, { 66.100000, 370.90000 }, { 67.300000, 373.30000 }, { 68.300000, 374.00000 }, { 70.000000, 374.40000 }, { 71.100000, 376.40000 }, { 72.900000, 377.80000 }, { 74.000000, 378.60000 }, { 75.000000, 379.60000 }, { 76.500000, 381.40000 }, { 77.000000, 425.40000 }, { 79.000000, 426.30000 }, { 80.000000, 426.80000 }, { 81.100000, 428.20000 }, { 82.300000, 429.30000 }, { 84.000000, 429.40000 }, { 85.100000, 430.30000 }, { 87.000000, 431.40000 } };


	if(setup_VUV_params) {
        fparslogNorm.SetParameters(landauNormpars);
	    fparsMPV.SetParameters(landauMPVpars);
	    fparsWidth.SetParameters(landauWidthpars);
	    fparsCte.SetParameters(expoCtepars);
	    fparsSlope.SetParameters(expoSlopepars);
	    //For far cases
	    fparslogNorm_far.SetParameters(landauNormpars_far);
	    fparsMPV_far.SetParameters(landauMPVpars_far);
	    fparsCte_far.SetParameters(expoCtepars_far);
	    //
	    setup_VUV_params = false;
	}

	std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
	//signals (remember this is transportation) no longer than 1us
	int i = (distance - 10.)/10.;
    
    double xmin = VUVTRange[i][0];
    double xmax = VUVTRange[i][1];
	const double signal_t_range = xmax;
	const double vuv_vgroup = 10.13;//cm/ns
	double t_direct = distance/vuv_vgroup;

	// Defining the two functions (Landau + Exponential) describing the timing vs distance
	double pars_landau[3];
	if(distance > d_break) {
		pars_landau[0]=fparsMPV_far.Eval(distance);
		pars_landau[1]=fparsWidth.Eval(d_break);
		pars_landau[2]=pow(10.,fparslogNorm_far.Eval(distance));
	}
	else {
	    pars_landau[0] = fparsMPV.Eval(distance);
	    pars_landau[1] = fparsWidth.Eval(distance);
	    pars_landau[2] = pow(10.,fparslogNorm.Eval(distance));
	}
	TF1 flandau ("flandau","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau.SetParameters(pars_landau);

	double pars_expo[2];
	if(distance > (d_break - 50.)) {
		pars_expo[0] = fparsCte_far.Eval(distance);
		pars_expo[1] = fparsSlope.Eval(d_break - 50.);
	}
	else {
	    pars_expo[0] = fparsCte.Eval(distance);
	    pars_expo[1] = fparsSlope.Eval(distance);
	}
	TF1 fexpo ("fexpo","expo",0, signal_t_range/2);
	fexpo.SetParameters(pars_expo);
    
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();
	//this is to find the intersection point between the two functions:
	TF1 fint ("fint",utility::finter_d,flandau.GetMaximumX(),3*t_direct,5);
	double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
	fint.SetParameters(parsInt);
	double t_int = fint.GetMinimumX();
	double minVal = fint.Eval(t_int);
	//the functions must intersect!!!
	if(minVal>0.015) {
		std::cout<<"WARNING: Parametrization of Direct Light discontinuous"<<std::endl;
	}
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();

	TF1 fVUVTiming ("fTiming",utility::LandauPlusExpoFinal,0,signal_t_range,6);
	double parsfinal[6] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
	fVUVTiming.SetParameters(parsfinal);
	
	// Set the number of points used to sample the function
	
	int f_sampling = 200;
	if(distance < 50)
		f_sampling *= 10;
	fVUVTiming.SetNpx(f_sampling);
	     
   	    
    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();
	for(int i=0; i<number_photons; i++) {
		double arrival_time = fVUVTiming.GetRandom(xmin,xmax);
		arrival_time_distrb.push_back(arrival_time);
	}

    
    std::chrono::steady_clock::time_point t6 = std::chrono::steady_clock::now();

    std::chrono::duration<double> time1 = std::chrono::duration_cast<std::chrono::duration<double>>(t1-t0);
    std::chrono::duration<double> time2 = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::chrono::duration<double> time3 = std::chrono::duration_cast<std::chrono::duration<double>>(t3-t2);
    std::chrono::duration<double> time4 = std::chrono::duration_cast<std::chrono::duration<double>>(t4-t3);
    std::chrono::duration<double> time5 = std::chrono::duration_cast<std::chrono::duration<double>>(t5-t4);
    std::chrono::duration<double> time6 = std::chrono::duration_cast<std::chrono::duration<double>>(t6-t5);
    
   // std::cout << "# = " << number_photons << "\t time1 = " << time1.count() << "\t time2 = " << time2.count() << "\t time3 = " << time3.count() << "\t time4 = " << time4.count() << "\t time5 = " << time5.count() << "\t time6 = " << time6.count() << "\t";
   // std::cout << "total = " << std::chrono::duration_cast<std::chrono::duration<double>>(t6-t0).count() << std::endl;
	//deleting ...
/*
        delete fparslogNorm;
        delete fparsMPV;
        delete fparsWidth;
        delete fparsCte;
        delete fparsSlope;
        delete flandau;
        delete fexpo;
        delete fint;
        delete fVUVTiming;

	return arrival_time_distrb;
}


//-----------------------------------------------------------
// Parametrization of the Visible light timing (result from direct transport + Rayleigh scattering ONLY)
// using a landau + expo function for a LDS configuration with TPC-coated foils ONLY covering the cathode.
// The function below returns the arrival time distribution given the distance IN CENTIMETERS between the
// scintillation/ionization point and the optical detector.
/*std::vector<double> utility::GetVisibleTimeOnlyCathode(double t0, int number_photons){
   std::vector<double> arrival_time_distrb;
   return arrival_time_distrb;
   }*/
/*
std::vector<double> utility::GetVisibleTimeOnlyCathode(double t0, int number_photons){
	//-----Distances in cm and times in ns-----//

	//gRandom->SetSeed(0);

	std::vector<double> arrival_time_distrb;
	arrival_time_distrb.clear();
	arrival_time_distrb.reserve(number_photons);
	// Parametrization data:
	double landauNormpars[4] = {7.54642, -0.441946, 0.0107579, -9.53399e-05};
	double landauMPVpars[4] = {-1.61482, 1.18624, 0.00105223, -9.52016e-05};
	double landauWidthpars[4] = {0.440124, -0.0557912, 0.00544957, -9.39128e-05};
	double expoCtepars[4] = {14.6874, -0.896761, 0.0214977, -0.000185728};
	double expoSlopepars[5] = {-0.650584, 0.0800897, -0.00379933, 7.91909e-05, -6.10836e-07};
	//range of t0s where the parametrization is valid [~8 - 55ns], then:
	const double t0_max = 55.;
	// Defining the function for the different parameters:
	TF1 fparslogNorm ("fparslogNorm","pol3",0,t0_max);
	fparslogNorm.SetParameters(landauNormpars);
	TF1 fparsMPV ("fparsMPV","pol3",0,t0_max);
	fparsMPV.SetParameters(landauMPVpars);
	TF1 fparsWidth ("fparsWidth","pol3",0,t0_max);
	fparsWidth.SetParameters(landauWidthpars);
	TF1 fparsCte ("fparsCte","pol3",0,t0_max);
	fparsCte.SetParameters(expoCtepars);
	TF1 fparsSlope ("fparsSlope","pol4",0,t0_max);
	fparsSlope.SetParameters(expoSlopepars);


	if(t0 < 8 || t0 > t0_max) {
		//std::cout<<"WARNING: Parametrization of Cathode-Only reflected Light not fully reliable"<<std::endl;
		//std::cout<<"Too close/far to the PMT  -> set 0 Visible photons(?)!!!!!!"<<std::endl;
		return arrival_time_distrb;
	}
	//signals (remember this is transportation) no longer than 1us
	const double signal_t_range = 1000.;

	double pars_landau[3] = {fparsMPV.Eval(t0), fparsWidth.Eval(t0),
		                 pow(10.,fparslogNorm.Eval(t0))};
	double pars_expo[2] = {fparsCte.Eval(t0), fparsSlope.Eval(t0)};
	double t0_break_point = 42.;//ns after the parametrization must be corrected (lack of statistics!)
	if(t0 > t0_break_point) {
		pars_landau[0] = -0.798934 + 1.06216*t0;
		pars_landau[1] = fparsWidth.Eval(t0_break_point);
		pars_landau[2] = pow(10.,fparslogNorm.Eval(t0_break_point));
		pars_expo[0] = fparsCte.Eval(t0_break_point);
		pars_expo[1] = fparsSlope.Eval(t0_break_point);
	}

	// Defining the two functions (Landau + Exponential) describing the timing vs t0
	TF1 flandau ("flandau","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau.SetParameters(pars_landau);
	TF1 fexpo ("fexpo","expo",0, signal_t_range/2);
	fexpo.SetParameters(pars_expo);
	//this is to find the intersection point between the two functions:
	TF1 fint ("fint",utility::finter_d,flandau.GetMaximumX(),2*t0,5);
	double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
	fint.SetParameters(parsInt);
	double t_int = fint.GetMinimumX();
	double minVal = fint.Eval(t_int);
	//the functions must intersect!!!
	//if(minVal>0.015)
	//std::cout<<"WARNING: Parametrization of Direct Light discontinuous (landau + expo)!!!!!!"<<std::endl;

	TF1 fVisTiming ("fTiming",utility::LandauPlusExpoFinal,0,signal_t_range,6);
	double parsfinal[6] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
	fVisTiming.SetParameters(parsfinal);
	// Set the number of points used to sample the function

	int f_sampling = 1000;
	if(t0 < 20)
		f_sampling *= 10;
	fVisTiming.SetNpx(f_sampling);

	for(int i=0; i<number_photons; i++)
		arrival_time_distrb.push_back(fVisTiming.GetRandom());

	//deleting ...
/*
        delete fparslogNorm;
        delete fparsMPV;
        delete fparsWidth;
        delete fparsCte;
        delete fparsSlope;
        delete flandau;
        delete fexpo;
        delete fint;
        delete fVisTiming;
 
	return arrival_time_distrb;
}

//--------------------------------------------------------------
// Parametrization of the Visible light timing (result from direct transport + Rayleigh scattering ONLY)
// using a landau + lndau function.The function below returns the arrival time distribution given the
// time of the first visible photon in the PMT, a weighted average of path lengths (in time) from the
// shortest reflection off each of the foil-coated walls, the distance IN CENTIMETERS between the scintillation/
// ionization point and the optical detectotr.
/*
   std::vector<double> utility::GetVisibleTimeFullConfig(double t0, double tmean, double distance, int number_photons) {
   std::vector<double> arrival_time_distrb;
   return arrival_time_distrb;
   }
 */
/*
std::vector<double> utility::GetVisibleTimeFullConfig2(double t0, double tmean, double distance, int number_photons) {

	std::vector<double> arrival_time_distrb;
	arrival_time_distrb.clear();
	arrival_time_distrb.reserve(number_photons);

	double divide = 0.43*tmean + 7;
	double vuv_vgroup = 10.13;//cm/ns
	double t_direct = distance/vuv_vgroup;

	std::string pols[6] = {"pol3", "pol1", "pol4", "pol2", "pol4", "pol2"};
	double landau1Norm[4] = {9.78924, -0.808646, 0.0286551, -0.000342326};
	double landau1MPV[2] = {-9.04501, 1.76972};
	double landau1Width[5] = {24.7515, -5.71531, 0.45703, -0.0144995, 0.000163086};
	double landau2Norm[3] = {3.44352, -0.0812814, 0.00118423};
	double landau2MPV[5] = {282.128, -57.8334, 4.50742, -0.143848, 0.00164436};
	double landau2Width[3] = {15.1667, -0.0786729, -0.000696796};

	///range of t0s where the parametrization is valid [~4 - 52ns (60ns for t_direct)], then:
	const double t0_max = 52.;
	const double t_direct_max = 60.;
	double t0_break_point = 35.;//ns after the parametrization must be corrected (lack of statistics!)
	double t0_max_mod = t_direct_max;

	if(t0 < 10 || t0 > t0_max || t_direct > t_direct_max) {
		//std::cout<<"WARNING: Parametrization of Full coverage reflected Light not fully reliable"<<std::endl;
		//std::cout<<"Too close/far to the PMT  -> set 0 Visible photons(?)!!!!!!"<<std::endl;
		return arrival_time_distrb;
	}


	// Defining the function for the different parameters:
	TF1 fparslogNorm1 ("fparslogNorm1",pols[0].c_str(),0,t0_max);
	fparslogNorm1.SetParameters(landau1Norm);
	TF1 fparsMPV1 ("fparsMPV1",pols[1].c_str(),0,t0_max);
	fparsMPV1.SetParameters(landau1MPV);
	TF1 fparsWidth1 ("fparsWidth1",pols[2].c_str(),0,t0_max);
	fparsWidth1.SetParameters(landau1Width);
	TF1 fparslogNorm2 ("fparslogNorm2",pols[3].c_str(),0,t0_max);
	fparslogNorm2.SetParameters(landau2Norm);
	TF1 fparsMPV2 ("fparsMPV2",pols[4].c_str(),0,t0_max);
	fparsMPV2.SetParameters(landau2MPV);
	TF1 fparsWidth2 ("fparsWidth2",pols[5].c_str(),0,t0_max_mod);
	fparsWidth2.SetParameters(landau2Width);


	//signals (remember this is transportation) no longer than 1us
	const double signal_t_range = 1000.;

	double pars_landau1[3] = {fparsMPV1.Eval(t0), fparsWidth1.Eval(t0), pow(10.,fparslogNorm1.Eval(t0))};
	double pars_landau2[3] = {fparsMPV2.Eval(t0), fparsWidth2.Eval(t_direct), pow(10.,fparslogNorm2.Eval(t0))};
	//To improve the degree of correlation in the estimation of Landau2-width,
	//in the population2 this parameter is estimated as a function of t_vuv, insted of t0

	if(t0 > t0_break_point) {
		pars_landau1[1] = fparsWidth1.Eval(t0_break_point);
		pars_landau1[2] = pow(10.,fparslogNorm1.Eval(t0_break_point));
		pars_landau2[0] = fparsMPV2.Eval(t0_break_point);
		pars_landau2[2] = pow(10.,fparslogNorm2.Eval(t0_break_point));
	}

	// Defining the two functions (Landau + Landau) describing the timing vs t0:
	TF1 flandau1 ("flandau1","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau1.SetParameters(pars_landau1);
	TF1 flandau2 ("flandau2","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau2.SetParameters(pars_landau2);
	TF1 fint ("fint",utility::finter_r,flandau1.GetMaximumX(),2*t0,6);
	double parsInt[6] = {pars_landau1[0], pars_landau1[1], pars_landau1[2], pars_landau2[0],
		             pars_landau2[1], pars_landau2[2]};
	fint.SetParameters(parsInt);
	double t_int = fint.GetMinimumX();
	double minVal = fint.Eval(t_int);

	TF1 fVisibleTiming ("fTiming",utility::LandauPlusLandauFinal,0,signal_t_range,7);
	double parsfinal[7] = {t_int, pars_landau1[0], pars_landau1[1], pars_landau1[2], pars_landau2[0],
		               pars_landau2[1], pars_landau2[2]};
	fVisibleTiming.SetParameters(parsfinal);
	// Set the number of points used to sample the function

	int f_sampling = 1500;
	fVisibleTiming.SetNpx(f_sampling);
	for(int i=0; i<number_photons; i++)
		arrival_time_distrb.push_back(fVisibleTiming.GetRandom());

	//deleting ...
/*
        delete fparslogNorm1;
        delete fparsMPV1;
        delete fparsWidth1;
        delete fparslogNorm2;
        delete fparsMPV2;
        delete fparsWidth2;
        delete flandau1;
        delete flandau2;
        delete fint;
        delete fVisibleTiming;

	return arrival_time_distrb;

}

std::vector<double> utility::GetVisibleTimeFullConfig1(double t0, double tmean, double distance, int number_photons) {
	//gRandom->SetSeed(0);

	std::vector<double> arrival_time_distrb;
	arrival_time_distrb.clear();
	arrival_time_distrb.reserve(number_photons);
	// Parametrization data:
	double divide = 0.43*tmean + 7;
	double vuv_vgroup = 10.13;//cm/ns
	double t_direct = distance/vuv_vgroup;

	// Defining the parameters of the two populations:
	std::string pols[6] = {"pol2", "pol2", "pol2", "pol2", "pol3", "pol2"};
	double landau1Norm[3] = {4.80632, -0.227272, 0.00409071};
	double landau1MPV[3] = {4.27391, 0.48747, 0.0312366};
	double landau1Width[3] = {0.789521, -0.0763977, 0.0094536};
	double landau2Norm[3] = {2.88774, -0.0188192, -0.00111117};
	double landau2MPV[4] = {-55.8751, 14.6612, -0.878218, 0.0198729};
	double landau2Width[3] = {10.5582, -0.539349, 0.0360326};

	///range of t0s where the parametrization is valid [~4 - 52ns (60ns for t_direct)], then:
	const double t0_max = 52.;
	const double t_direct_max = 60.;
	double t0_break_point = 35.;//ns after the parametrization must be corrected (lack of statistics!)
	double t0_max_mod = t0_max;

	if(t0 < 4 || t0 > t0_max) {
		//std::cout<<"WARNING: Parametrization of Full coverage reflected Light not fully reliable"<<std::endl;
		//std::cout<<"Too close/far to the PMT  -> set 0 Visible photons(?)!!!!!!"<<std::endl;
		return arrival_time_distrb;
	}

	// Defining the function for the different parameters:
	TF1 fparslogNorm1 ("fparslogNorm1",pols[0].c_str(),0,t0_max);
	fparslogNorm1.SetParameters(landau1Norm);
	TF1 fparsMPV1 ("fparsMPV1",pols[1].c_str(),0,t0_max);
	fparsMPV1.SetParameters(landau1MPV);
	TF1 fparsWidth1 ("fparsWidth1",pols[2].c_str(),0,t0_max);
	fparsWidth1.SetParameters(landau1Width);
	TF1 fparslogNorm2 ("fparslogNorm2",pols[3].c_str(),0,t0_max);
	fparslogNorm2.SetParameters(landau2Norm);
	TF1 fparsMPV2 ("fparsMPV2",pols[4].c_str(),0,t0_max);
	fparsMPV2.SetParameters(landau2MPV);
	TF1 fparsWidth2 ("fparsWidth2",pols[5].c_str(),0,t0_max_mod);
	fparsWidth2.SetParameters(landau2Width);


	//signals (remember this is transportation) no longer than 1us
	const double signal_t_range = 1000.;

	double pars_landau1[3] = {fparsMPV1.Eval(t0), fparsWidth1.Eval(t0), pow(10.,fparslogNorm1.Eval(t0))};
	double pars_landau2[3] = {fparsMPV2.Eval(t0), fparsWidth2.Eval(t0), pow(10.,fparslogNorm2.Eval(t0))};
	//To improve the degree of correlation in the estimation of Landau2-width,
	//in the population2 this parameter is estimated as a function of t_vuv, insted of t0
	//prevent "bad" behaviour (below 6 ns) of the pol parametrizing this parameter (lack in stats):
	if(t0 < 6)
		pars_landau2[0] = fparsMPV2.Eval(6.);

	// Defining the two functions (Landau + Landau) describing the timing vs t0:
	TF1 flandau1 ("flandau1","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau1.SetParameters(pars_landau1);
	TF1 flandau2 ("flandau2","[2]*TMath::Landau(x,[0],[1])",0,signal_t_range/2);
	flandau2.SetParameters(pars_landau2);
	TF1 fint ("fint",utility::finter_r,flandau1.GetMaximumX(),2*t0,6);
	double parsInt[6] = {pars_landau1[0], pars_landau1[1], pars_landau1[2], pars_landau2[0],
		             pars_landau2[1], pars_landau2[2]};
	fint.SetParameters(parsInt);
	double t_int = fint.GetMinimumX();
	double minVal = fint.Eval(t_int);
	//the functions must intersect!!!
	//if(minVal>0.015)
	//std::cout<<"WARNING: Parametrization of Reflected Light discontinuous (landau + landau)!!!!!!"<<std::endl;

	TF1 fVisibleTiming ("fTiming",utility::LandauPlusLandauFinal,0,signal_t_range,7);
	double parsfinal[7] = {t_int, pars_landau1[0], pars_landau1[1], pars_landau1[2], pars_landau2[0],
		               pars_landau2[1], pars_landau2[2]};
	fVisibleTiming.SetParameters(parsfinal);
	// Set the number of points used to sample the function

	int f_sampling = 1500;
	fVisibleTiming.SetNpx(f_sampling);
	for(int i=0; i<number_photons; i++)
		arrival_time_distrb.push_back(fVisibleTiming.GetRandom());

	//deleting ...
/*
        delete fparslogNorm1;
        delete fparsMPV1;
        delete fparsWidth1;
        delete fparslogNorm2;
        delete fparsMPV2;
        delete fparsWidth2;
        delete flandau1;
        delete flandau2;
        delete fint;
        delete fVisibleTiming;

	return arrival_time_distrb;
}

//----------------------------------------------------
// Get random time (ns) from a parametrized approximation to the reflected
// light distribution, as obtained in a toy MC.
// (UNDER CONSTRUCTION)
// Input argument units: meters
// WARNING: Used to get the the "weighted mean time" to be used by the 2-Landau reflected
// light parametrization, not the time distribution!
// ---------------
double utility::TimingParamReflected(TVector3 ScintPoint, TVector3 OpDetPoint ) {

	// This function takes as input the scintillation point and the optical
	// detector point (in m), and uses this information to return a photon propagation
	// time assuming the field cage walls and cathode are covered ~100% uniformly
	// in TPB reflector foil (100% efficiency for VUV conversion, 95% diffuse
	// reflectance to blue-visible light).
	//
	// Note:  Modeling photon arrival times as a function of scintillation
	//        position in a fully reflective system is complicated, given there
	//        is no obvious "distance to PMT" handle we can use to parametrize
	//        the distribution.  This parametrization is only an approximation
	//        to this timing, and should be treated accordingly.
	//
	//        The arrival time distributions are well-described by the double-
	//        landau function, DoubleLandauCutoff.  However, this has several
	//        free parameters that were not able to be parametrized.  So for
	//        now, we are zeroing-out the second ("fast") landau in the function
	//        and thus treating it as a single landau with a t0-cutoff.
	//
	//        This parametrization was developed specifically for the SBND
	//        half-TPC (2x4x5m).  It's not yet known if it holds for other
	//        dimensions.  These dimensions are thus hard-coded in for now.
	//

	// Convert to center-of-detector origin coordinates, reversing the
	// x direction to go from APA to cathode (RHR!!)
	ScintPoint[0] = -1.*(ScintPoint[0]-1.);
	OpDetPoint[0] = -1.*(OpDetPoint[0]-1.);
	ScintPoint[2] -= 2.5;
	OpDetPoint[2] -= 2.5;

	// Each plane's position in x,y,z
	//    index 0: APA (x)
	//    index 1: cathode plane (x)
	//    index 2: bottom (y)
	//    index 3: top (y)
	//    index 4: upstream wall (z)
	//    index 5: downstream wall (z)
	//double  plane_depth[6] = {-1.,1.,-2.,2.,-2.5,2.5};
	// to match the position of the foils in the gdml file
	double plane_depth[6] = {-1.,1.,-2.-0.015,2.+0.015,-2.5-0.015,2.5+0.015};
	int dir_index_norm[6] = {0,0,1,1,2,2};

	// Speed of light in LAr, refractive indices
	double n_LAr_VUV = 2.632;       // Effective index due to group vel.
	double n_LAr_vis = 1.23;
	double c_LAr_VUV = 0.12;        // m per s
	double c_LAr_vis = 0.29979/n_LAr_vis; // m per s

	// TVectors to be used later
	TVector3 image(0,0,0);
	TVector3 bounce_point(0,0,0);
	TVector3 hotspot(0,0,0);
	TVector3 v_to_wall(0,0,0);

	// First find shortest 1-bounce path and fill
	// vector of "tAB" paths, each with their weighting:
	double t0   = 99;
	double tAB[5];
	double WAB=0;
	double tAB_sum=0;
	int counter = 0;
	double d0 = 1000;
	for(int j = 1; j<6; j++) {

		v_to_wall[dir_index_norm[j]] = plane_depth[j] - ScintPoint[dir_index_norm[j]];

		// hotspot is point on wall where TPB is
		// activated most intensely by the scintillation
		hotspot = ScintPoint + v_to_wall;

		// define "image" by reflecting over plane
		image = hotspot + v_to_wall*(n_LAr_vis/n_LAr_VUV);

		// find point of intersection with plane j of ray
		// from the PMT to the image
		TVector3 tempvec = (OpDetPoint-image).Unit();
		double tempnorm= ((image-hotspot).Mag())/tempvec[dir_index_norm[j]];
		bounce_point = image + tempvec*tempnorm;

		// find t, check for t0
		double t1 = (ScintPoint-bounce_point).Mag()/c_LAr_VUV;
		double t2 = (OpDetPoint-bounce_point).Mag()/c_LAr_vis;
		double t  = t1 + t2;
		double d = (ScintPoint-bounce_point).Mag() + (OpDetPoint-bounce_point).Mag();
		if( t < t0 ) t0 = t;
		if( d < d0 ) d0 = d;

		// find "tAB" and its weight
		double dA = (ScintPoint-hotspot).Mag();
		double dB = (OpDetPoint-hotspot).Mag();
		double tA = dA/c_LAr_VUV;
		double tB = dB/c_LAr_vis;
		double hotspot_weight = 1./pow(dA,2.) - 0.0294/pow(dA,3.);
		double tAB_w          = hotspot_weight/(1.+dB*dB);
		tAB[counter]  =  tA + tB;
		tAB_sum       += (tA+tB)*tAB_w;
		WAB           += tAB_w;

		counter++;

	} //<-- end loop over 5 foil-covered planes

	// weighted mean
	double tAB_mean = tAB_sum/WAB;

	// now find tAB_spread
	double tAB_spread = 0;
	for(int jj=0; jj<counter; jj++) tAB_spread += fabs( tAB[jj] - tAB_mean )/counter;

	return tAB_mean;
	//static double t_Will[2] = {tAB_mean, t0};

	//return t_Will;


	/*
	   // Now we have all the info we need for the parametrization
	   double mpv    = 93.695 - 16.572*tAB_spread + 2.064*pow(tAB_spread,2) - 0.1177*pow(tAB_spread,3) + 0.002025*pow(tAB_spread,4);
	   double width  = 14.75 - 0.4295*tAB_spread;

	   TF1 *fReflectedTiming = new TF1("fReflectedTiming",DoubleLandauCutoff,0.,250.,7);
	   fReflectedTiming -> SetNpx(10000);
	   fReflectedTiming -> FixParameter(0,t0);
	   fReflectedTiming -> FixParameter(1,mpv);
	   fReflectedTiming -> FixParameter(2,width);
	   fReflectedTiming -> FixParameter(3,0);  // ***setting "fast" landau to zero***
	   fReflectedTiming -> FixParameter(4,1.); // fast landau mpv (irrelevant for now)
	   fReflectedTiming -> FixParameter(5,1.); // fast landau width (irrelevant for now)
	   fReflectedTiming -> FixParameter(6,1.); // overall normalization

	   return fReflectedTiming->GetRandom();
	 
}

std::vector<double> utility::TimingParamReflected2(TVector3 ScintPoint, TVector3 OpDetPoint ) {
    //// ALTERED FOR DUNE: Removed multiple reflective walls, as only cathode is
    //// is relevant here.
    //
	// This function takes as input the scintillation point and the optical
	// detector point (in m), and uses this information to return a photon propagation
	// time assuming the field cage walls and cathode are covered ~100% uniformly
	// in TPB reflector foil (100% efficiency for VUV conversion, 95% diffuse
	// reflectance to blue-visible light).
	//
	// Note:  Modeling photon arrival times as a function of scintillation
	//        position in a fully reflective system is complicated, given there
	//        is no obvious "distance to PMT" handle we can use to parametrize
	//        the distribution.  This parametrization is only an approximation
	//        to this timing, and should be treated accordingly.
	//
	//        The arrival time distributions are well-described by the double-
	//        landau function, DoubleLandauCutoff.  However, this has several
	//        free parameters that were not able to be parametrized.  So for
	//        now, we are zeroing-out the second ("fast") landau in the function
	//        and thus treating it as a single landau with a t0-cutoff.
	//
	//        This parametrization was developed specifically for the SBND
	//        half-TPC (2x4x5m).  It's not yet known if it holds for other
	//        dimensions.  These dimensions are thus hard-coded in for now.
	//

	// to match the position of the foils in the gdml file
	double plane_depth = 363.38405; //NOT CONFIRMED YET

	// Speed of light in LAr, refractive indices
	double n_LAr_VUV = 2.632;       // Effective index due to group vel.
	double n_LAr_vis = 1.23;
	double c_LAr_VUV = 0.12;        // m per s
	double c_LAr_vis = 0.29979/n_LAr_vis; // m per s

	// TVectors to be used later
	TVector3 image(0,0,0);
	TVector3 bounce_point(0,0,0);
	TVector3 hotspot(0,0,0);
	TVector3 v_to_wall(0,0,0);

	// First find shortest 1-bounce path and fill
	// vector of "tAB" paths, each with their weighting:
	//double t0   = 99;
	double tAB;
	//double d0 = 1000;
	
	v_to_wall[0] = plane_depth - ScintPoint[0];

	// hotspot is point on wall where TPB is
	// activated most intensely by the scintillation
	hotspot = ScintPoint + v_to_wall;

	// define "image" by reflecting over plane
	image = hotspot + v_to_wall*(n_LAr_vis/n_LAr_VUV);

	// find point of intersection with plane j of ray
	// from the PMT to the image
	TVector3 tempvec = (OpDetPoint-image).Unit();
	double tempnorm= ((image-hotspot).Mag())/abs(tempvec[0]); //added abs() to denom. to ensure +ve norm
	bounce_point = image + tempvec*tempnorm;

	// find t, check for t0
	double t1 = (ScintPoint-bounce_point).Mag()/c_LAr_VUV;
	double t2 = (OpDetPoint-bounce_point).Mag()/c_LAr_vis;
	double t  = t1 + t2;
	double d = (ScintPoint-bounce_point).Mag() + (OpDetPoint-bounce_point).Mag();

	double dA = (ScintPoint-hotspot).Mag();
	double dB = (OpDetPoint-hotspot).Mag();
	double tA = dA/c_LAr_VUV;
	double tB = dB/c_LAr_vis;
	tAB  =  tA + tB;
	
	std::vector<double> t_Will= {tAB, t};
    
    /*cout << "ScintPoint:\t";        ScintPoint.Print();
    cout << "OpDetPoint:\t";        OpDetPoint.Print();
    cout << "image:\t\t";           image.Print();
    cout << "tempvec:\t";           tempvec.Print();
    cout << "bounce_point:\t";      bounce_point.Print();
    cout << "hotspot:\t";           hotspot.Print();
	return t_Will;
}

TVector3 utility::GetShortestPathPoint(TVector3 ScintPoint, TVector3 OpDetPoint) {
    /*
    Determines the position of X, the point at which at photon is reflected on a shortest-time path.
    For the DUNE detector exclusively.
    
    // to match the position of the foils in the gdml file

	double plane_depth = 363.38405; //NOT CONFIRMED YET

	// Speed of light in LAr, refractive indices
	double n_LAr_VUV = 2.632;       // effective index due to group vel.
	double n_LAr_vis = 1.23;
	double c_LAr_VUV = 0.12;        // m per s
	double c_LAr_vis = 0.29979/n_LAr_vis; // m per s

	// TVectors to be used later
	TVector3 image(0,0,0);
	TVector3 bounce_point(0,0,0);
	TVector3 hotspot(0,0,0);
	TVector3 v_to_wall(0,0,0);

	// First find shortest 1-bounce path and fill
	// vector of "tAB" paths, each with their weighting:
	//double t0   = 99;
	double tAB;
	//double d0 = 1000;
	
	v_to_wall[0] = plane_depth - ScintPoint[0];
	
	// hotspot is point on wall where TPB is
	// activated most intensely by the scintillation
	hotspot = ScintPoint + v_to_wall;

	// define "image" by reflecting over plane
	image = hotspot + v_to_wall*(n_LAr_vis/n_LAr_VUV);

	// find point of intersection with plane j of ray
	// from the PMT to the image
	TVector3 tempvec = (OpDetPoint-image).Unit();
	double tempnorm= ((image-hotspot).Mag())/std::abs(tempvec[0]); //added abs() to denom. to ensure +ve norm
	bounce_point = image + tempvec*tempnorm;
    
    return bounce_point;
}

TVector3 utility::GetShortestPathPoint_2(TVector3 ScintPoint, TVector3 OpDetPoint){
	// determines and returns position of point at which photon is reflected on the fastest time path
	// method 2: minimising VUV path length, since Rayleigh scattering substantial for this part but negligible for refl light path

	// plane depth, mathcing position of foils in gdml file
	double plane_depth = 363.38405;

	// calculate bounce point, shortest path just straight line to foils position
	TVector3 bounce_point(plane_depth,ScintPoint[1],ScintPoint[2]);
	
	return bounce_point; 
}

std::vector<double> utility::GetReflTime(double distance, int number_photons) {
    /*
    Outputs a 1D vector array of times for each produced photon on the visible part of its journey
    
    std::vector<double> ReflTimes(number_photons,0);
    double mean_v = 23.99; // cm/ns
    double rms_v = 0.07508; // cm/ns
    double v;
    
    for (int i=0; i<number_photons; i++) {
        v = gRandom->Gaus(mean_v,rms_v);
        ReflTimes[i] = distance/v;
    }
    
    return ReflTimes;
}

std::vector<double> utility::GetVisTime0thOrder(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons) {
    /*
    Attempts to perform a zeroth-order parametrization of the time distribution for reflected light in DUNE.
    Looks only at the shortest-time path for light, where there are two contributions:
     - VUV light towards the cathode; time through GetVUVTime() function
     - reflected visible light back towards APA; distance constant but velocity tightly Gaussian-distributed
    INPUTS:
     - TVector3 ScintPoint: position of scintillation in DUNE (in cm)
     - TVector3 OpDetPoint: position of photon detector in DUNE (in cm)
     - number_photons:      number of photons produced in scintillation
    OUTPUT:
     - 1D vector array of times (in ns) for each produced photon
    
   
   	// get point of reflection for shortest path
    TVector3 X = utility::GetShortestPathPoint(ScintPoint, OpDetPoint);
    
    // calculate distances travelled by each part
    double VUVdist = (X-ScintPoint).Mag();
    double Visdist = (OpDetPoint-X).Mag();
    
    std::vector<double> VUVTimes  = utility::GetVUVTime (VUVdist, number_photons);
    std::vector<double> ReflTimes = utility::GetReflTime(Visdist, number_photons);

    std::vector<double> TotTimes(number_photons,0);
    for (int i=0; i<number_photons; i++) {
        TotTimes[i] = VUVTimes[i] + ReflTimes[i];
        //std::cout << "VUVtime " << i << " = " << VUVTimes[i] << "	ReflTime " << i << " = " << ReflTimes[i] << std::endl;
    }
    return TotTimes;
}

*/