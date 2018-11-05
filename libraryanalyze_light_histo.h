#ifndef LIBRARYANALYZE_LIGHT_HISTO_H
#define LIBRARYANALYZE_LIGHT_HISTO_H

#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "THStack.h"
#include "TColor.h"
#include "TLegend.h"
#include "TMarker.h"

// -----------------------------------------------------------
// -----------------------------------------------------------
// Use this file to set the configurable parameters for the mc
// -----------------------------------------------------------
// -----------------------------------------------------------

using namespace std;

///-------------------------------------
/// Custom Sorting Function
///-------------------------------------
bool sort_function(std::pair<double, int> pair1, std::pair<double, int> pair2)
{ return (pair1.first < pair2.first); }
///-------------------------------------


///-------------------------------------
//--------WHAT to generate?-------------
///-------------------------------------
bool fixed_energy = true; double fixedE = 20.0; //MeV
bool supernova = false;
bool gen_argon = false;
bool gen_radon = false;
///-------------------------------------
///-------------------------------------
///-------DO timing calculations?-------
///-------------------------------------
bool do_timings = true;
double step_size = 1.0; 	// step size for discretisation of timing array in cm
///-------------------------------------
//--------WHERE to generate?-------------
///-------------------------------------
// Choose one only!
bool random_pos = true;	// works
double PosMin[3] = {10,-600,300}; 	//For random_pos option, generate in this range
double PosMax[3] = {330,600,1000};
bool fixed_xpos = false; 	// needs updating, range getting random position from is not valid for dune library
bool fixed_pos = false;		// works
double fixedX = 250; 		// cm 
double fixedY = 31.1784; 	// cm 
double fixedZ = 580.099; 	// cm


///-------------------------------------
//-------- make time cut?-------------
///-------------------------------------
bool cut = false; 		// NB you can always make time cuts when you're analysing the files
double time_cut = 0.1; 	// in microseconds - 0.1 mu_s = 100 ns


///-------------------------------------
//------Light System Configuration------
///-------------------------------------
//--------------------------------------
///Config List:
///0 = Full Foils
///1 = Cath Foils
///2 = VUV only
const int config = 1; // cathode foils only configuration for dune
//--------------------------------------
//--------------------------------------

//--------------------------------------
// These bools are used in the scintillation timing functions
std::string libraryfile;
bool reflected;
bool reflT;
//--------------------------------------

//--------------------------------------
//TTree branches and data products:
//-------------------------------------
TFile event_file("event_file_testing.root", "RECREATE", "Event File");

TTree *data_tree = new TTree("data_tree", "data tree");
TTree *data_tree_vuv = new TTree("data_tree_vuv", "data tree_vuv");
TTree *data_tree_vis = new TTree("data_tree_vis", "data tree_vis");

TTree *event_tree = new TTree("event_tree", "event tree");
double data_time;
double data_time_vuv;
double data_time_vis;

int data_pmt;
int data_pmt_vuv;
int data_pmt_vis;

int data_event;
int data_event_vuv;
int data_event_vis;

double data_x_pos;
double data_x_pos_vuv;
double data_x_pos_vis;

double data_y_pos;
double data_y_pos_vuv;
double data_y_pos_vis;

double data_z_pos;
double data_z_pos_vuv;
double data_z_pos_vis;

int event_no;
int event_vox;
double event_x_pos;
double event_y_pos;
double event_z_pos;
double event_E;
//--------------------------------------

///-------------------------------------
//----LAr & Ar39 properties--------------------
///-------------------------------------
const double MassE = 0.510998910; 	// mass electron - MeV/c^2
const double Q_Ar = 0.565;			//Q value of decay - Ar39

const double t_singlet = 0.000000006; 		//6ns
const double t_triplet = 0.0000015; 		//1.5 us
const double scint_time_window = 0.00001; 	//10 us

const int scint_yield_electron = 24000;	// Scintillation yield of LAr at 500 V/cm
const double activity_Ar = 1.; 			// Ar39 roughly 1 Bq/k
///-------------------------------------

///-------------------------------------
//----SN properties---------------------
///-------------------------------------
const double Eav = 20.;				// Average energy for SN spectrum
const double expected_sn = 2.8;		// For poisson weighting
///-------------------------------------

///-------------------------------------
//----Radon properties---------------------
///-------------------------------------
const int scint_yield_alpha = 16800; 	// SY of alpha particles at 500 V/cm - from larsoft
const double activity_Rn = 0.000021; 	// Bq/kg
const double massAlpha = 3727.3794; 	// alpha particle mass - MeV/c^2
const double Q_Rn = 5.590; 				// deposited energy from a radon decay - Rn-222 --> Po-218
///-------------------------------------

///-------------------------------------
//----TPC and PMT properties---------------------
///-------------------------------------
// QE of Arapucas, mesh efficiency  + /0.46 to remove bar attenuation factor included in library
const double quantum_efficiency = 0.01*0.7; //0.46; 		
const double catcov = 0.8; 	// proportion of cathode covered by TPB
const double vuvfrac = 0.4;
const double visfrac = 1.; 	// inclusive mode for now
const double mass = 112000.; //SBND 112ton LAr
const double time_window = 10.; //(0.0012 * 10.);//1.2 [ms] is the readout window.
const double time_frames = time_window/0.0012;
///-------------------------------------

///-------------------------------------
//----Number of Events------------------
///-------------------------------------
// Fixed energy (electron like) events:
const int max_events_FE = 100000;

// Ar-39 events:
//const int max_events_Ar = 10;
const int max_events_Ar = activity_Ar * mass/2 * time_window;//Half volume for 1 TPC
const int Ar_decays_per_sec = activity_Ar* mass/2; // decay rate in one TPC

// Radon events:
const int max_events_Rn = 1;
//const int max_events_Rn = activity_Rn * mass/2 * time_window;//Half volume for 1 (NOTE: for a small time window, this will probably return 0)
const double Rn_decays_per_sec = activity_Rn* mass/2; // decay rate in one TPC

// Supernova events:
const int max_events_SN = 1;
//int max_events_SN = utility::poisson(expected_sn,gRandom->Uniform(1.),1.);


//--------------------------------------
//--------------------------------------
// Don't change these
// used to read optical channel locations into
vector<vector<double>> myfile_data;
// initialise library access class
LibraryAccess lar_light;
// lists of variables for generating
vector<double> energy_list;
vector<double> decay_time_list;
vector<vector<double>> position_list;
vector<int> voxel_list;

#endif
