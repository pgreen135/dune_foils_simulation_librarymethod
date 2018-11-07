#ifndef TIMEPARAM_H
#define TIMEPARAM_H

// class containing imformation required for VUV timings parameterisation
// builds array of parameterisations for pre-determined distances
// calculates and returns the VUV transport time distribution using pre-built array of parameterisations

// also contains function calculating vis light transport times
// calculates shortest path considering refractive indices of each component in LAr
// uses smearing of VUV parameterisation distribution to approximate the vis transport time distribution

#include <vector>

#include "TF1.h"
#include "TVector3.h"


class timeparam {

private:
	// *************************************************************************************************
	// Discretisation of VUV time parameterisations
	// *************************************************************************************************
	// discretisation step size in cm, set in libraryanalyze_light_histo.h
	double step_size;
	// maximum distance in cm parameterisations are generated for (only generated when required)
	double d_max = 2500;
	// vector containing generated VUV timing parameterisations
	std::vector<TF1> VUV_timing;
	// vector containing min and max range VUV timing parameterisations are sampled to
	std::vector<double> VUV_max;
	std::vector<double> VUV_min;
	
	// *************************************************************************************************
	// Definition of the parameters for VUV timings parameterisation
	// *************************************************************************************************
	// Parameters of the Landau + Exponential (<= 300 cm) and Landau (> 300 cm) models
	// Landau parameters 
	std::vector<double> vDistances_all = {62.5, 87.5, 112.5, 137.5, 162.5, 187.5, 212.5, 237.5, 262.5, 287.5, 312.5, 337.5, 362.5, 387.5, 412.5, 437.5, 462.5, 487.5, 512.5, 537.5, 562.5, 587.5, 612.5, 637.5, 662.5, 687.5, 712.5, 737.5, 762.5, 787.5, 812.5, 837.5, 862.5, 887.5, 912.5, 937.5, 962.5, 987.5, 1012.5};
	std::vector<double> vNorm_over_entries = {1.46739, 0.880294, 0.571107, 0.344403, 0.261388, 0.20993, 0.170956, 0.152152, 0.136094, 0.120939, 0.109945, 0.10561, 0.0958589, 0.0904341, 0.0839861, 0.0803693, 0.0761224, 0.0716448, 0.0689223, 0.0668696, 0.0635322, 0.0616896, 0.0600554, 0.0577425, 0.0561707, 0.0539962, 0.0517729, 0.050872, 0.0488119, 0.0469635, 0.0466058, 0.0453283, 0.04366, 0.0425648, 0.0411353, 0.0396493, 0.0388351, 0.0383864, 0.0385728};
	std::vector<double> vMpv = {6.05668, 9.05062, 12.3438, 17.1553, 21.9225, 26.9042, 33.1595, 37.1543, 43.6424, 49.9016, 55.7976, 59.7093, 66.7296, 72.5434, 78.6925, 84.1388, 90.0175, 97.1053, 101.911, 107.799, 113.748, 119.527, 124.167, 130.438, 134.988, 141.886, 147.848, 153.056, 160.653, 168.234, 171.821, 178.68, 186.341, 193.222, 202.61, 212.357, 222.014, 232.012, 235.169};
	std::vector<double> vWidth = {0.977654, 1.57296, 2.42236, 4.15371, 5.76189, 7.49866, 9.87714, 10.9349, 13.4434, 15.5028, 19.8005, 20.6927, 22.9557, 24.5475, 26.2879, 27.7879, 29.296, 31.1889, 32.5755, 33.5919, 35.4932, 36.7544, 37.8536, 39.5658, 40.8131, 42.6414, 44.6617, 45.9231, 48.2682, 50.799, 51.9537, 54.1583, 57.8934, 61.4115, 65.4056, 70.7218, 74.6922, 79.6439, 81.8717};
	// Exponential parameters
	std::vector<double> vDistances = {55, 65, 75, 85, 95, 105, 115, 125, 135, 145, 155, 165, 175, 185, 195, 205, 215, 225, 235, 245, 255, 265, 275, 285, 295};
	std::vector<double> vSlope = {-0.0885849, -0.0688936, -0.0519349, -0.0481695, -0.0408145, -0.0359584, -0.0341716, -0.0302902, -0.0282722, -0.026348, -0.0240983, -0.023172, -0.0220861, -0.0207058, -0.0197934, -0.0193794, -0.0182472, -0.0175204, -0.0173492, -0.0163832, -0.0164717, -0.0160823,-0.0153199, -0.0149124, -0.0142298};
	// Line fit to the profiles of the ratio between the normalization parameters of the exponential and the landau functions vs distance (<= 300 cm)
	// Fits made to profiles in three different offset angle bins [0, 30deg], [30deg, 60deg] and [60deg, 90 deg]
	const double Expo_over_Landau_norm[3][2] = { {-4.44152e-02, 1.25321e-03}, {-2.74406e-02, 1.27240e-03}, {3.48944e-02, 1.27080e-03} };
	// VUV group velocity
	const double vuv_vgroup_mean = 10.13;//cm/ns
	const double vuv_vgroup_max = 15.;//cm/ns
	// Distance from "Landau + Expo" -> "Single Landau" model
	const double inflexion_point_distance = 300.; //cm

	// *************************************************************************************************
	// Definition of the parameters for vis timings calculation
	// *************************************************************************************************
	// velocity of visible light
	const double vis_vmean = 23.99; // cm/ns
	const double vis_vrms = 0.07508; // cm/ns
	// x-position of cathode foils, matching position of foils in gdml file - this needs changing if using different geometry
	const double plane_depth = 363.38405;
	// refractive indices in LAr
	double n_LAr_VUV = 2.632;       // effective index due to group vel.
	double n_LAr_vis = 1.23;
	double c_LAr_VUV = 0.12;        		// m/s
	double c_LAr_vis = 0.29979/n_LAr_vis; 	// m/s

	double pi =  3.141592653589793;


public:
	// constructor
	timeparam(double step);

	// destructor
	~timeparam(){std::cout << "Deleting timing parameterisations.\n";};

	// parameterisation generation function
	void generateparam(int index);

    // VUV arrival times calculation function
    std::vector<double> getVUVTime(double distance, int number_photons);

    // vis arrival times calculation function
    std::vector<double> getVisTime(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons);
};

#endif