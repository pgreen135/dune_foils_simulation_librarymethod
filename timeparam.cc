// implementation of timeparam class

#include "timeparam.h"
#include "timeparamservices.h"

#include <string>
#include <cmath>

#include "TMath.h"
#include "TRandom.h"

using namespace std;

// constructor
timeparam::timeparam(double size): step_size{size} {	
	
    // create empty parameterisations array, parameterisations generated as they are required
    int num_params = (d_max - 25) / step_size;  // for d < 25cm, no parameterisaton - delta function instead
	// default TF1() constructor gives function with 0 dimensions, can then check numDim to see if parameterisation has been generated	
	vector<TF1> VUV_timing_temp(num_params,TF1());
	VUV_timing = VUV_timing_temp;
    // initialise vectors to contain range parameterisations sampled to in each case
    vector<double> VUV_empty(num_params, 0);
    VUV_max = VUV_empty;
    VUV_min = VUV_empty;
}

// parameterisation generation function
void timeparam::generateparam(int index) {
    gRandom->SetSeed(0);

    // get distance 
    double distance_in_cm = (index * step_size) + 25;
    
    // time range
    const double signal_t_range = 5000.;
   
    // parameterisation TF1    
    TF1* fVUVTiming;
      
    // For very short distances the time correction is just a shift
    double t_direct_mean = distance_in_cm/vuv_vgroup_mean;
    double t_direct_min = distance_in_cm/vuv_vgroup_max;
      
    // Defining the model function(s) describing the photon transportation timing vs distance 
    // Getting the landau parameters from the time parametrization
    double* pars_landau = interpolate(vDistances_all, vMpv, vWidth, vNorm_over_entries, distance_in_cm, true);
    // Deciding which time model to use (depends on the distance)
    // defining useful times for the VUV arrival time shapes
    if(distance_in_cm >= inflexion_point_distance) {
        double pars_far[4] = {t_direct_min, pars_landau[0], pars_landau[1], pars_landau[2]};
        // Set model: Landau 
        fVUVTiming =  new TF1("fVUVTiming",model_far,0,signal_t_range,4);
        fVUVTiming->SetParameters(pars_far);
    }
    else {
        // Set model: Landau + Exponential 
        fVUVTiming =  new TF1("fVUVTiming",model_close,0,signal_t_range,7); 
        // Exponential parameters
        double pars_expo[2];   
        // Getting the exponential parameters from the time parametrization
        pars_expo[1] = interpolate(vDistances, vSlope, distance_in_cm, true);
        pars_expo[0] = Expo_over_Landau_norm[1][0] + Expo_over_Landau_norm[1][1]*distance_in_cm;
        pars_expo[0] *= pars_landau[2];
        pars_expo[0] = log(pars_expo[0]);
        // this is to find the intersection point between the two functions:
        TF1* fint = new TF1("fint",finter_d,pars_landau[0],4*t_direct_mean,5);
        double parsInt[5] = {pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1]};
        fint->SetParameters(parsInt);
        double t_int = fint->GetMinimumX();
        double minVal = fint->Eval(t_int);
        // the functions must intersect - output warning if they don't
        if(minVal>0.015) {
            std::cout<<"WARNING: Parametrization of VUV light discontinuous for distance = " << distance_in_cm << std::endl;
        }
        delete fint;   
        double parsfinal[7] = {t_int, pars_landau[0], pars_landau[1], pars_landau[2], pars_expo[0], pars_expo[1], t_direct_min};
        fVUVTiming->SetParameters(parsfinal);    
    }

    // calculate max and min distance relevant to sample parameterisation 
    // max
    const int nq_max=1;
    double xq_max[nq_max];
    double yq_max[nq_max];    
    xq_max[0] = 0.99;
    fVUVTiming->GetQuantiles(nq_max,yq_max,xq_max);
    double max = yq_max[0];
    // min
    double min = t_direct_min;

    // set the number of points used to sample parameterisation
    // for shorter distances, peak is sharper so more sensitive sampling required - values could be optimised 
    int f_sampling;
    if (distance_in_cm < 50) { f_sampling = 10000; }
    else if (distance_in_cm < 100){ f_sampling = 5000; }
    else { f_sampling = 1000; }
    fVUVTiming->SetNpx(f_sampling);    

    // generate the sampling
    // the first call of GetRandom generates the timing sampling and stores it in the TF1 object, this is the slow part
    // all subsequent calls check if it has been generated previously and are quick
    double arrival_time = fVUVTiming->GetRandom(min,max);
    // add timing to the vector of timings and range to vectors of ranges
    VUV_timing[index] = *fVUVTiming;
    VUV_max[index] = max;
    VUV_min[index] = min;

    delete fVUVTiming;
}

// VUV arrival times calculation function
vector<double> timeparam::getVUVTime(double distance, int number_photons) {
    gRandom->SetSeed(0);

    // pre-allocate memory
    std::vector<double> arrival_time_distrb;
    arrival_time_distrb.clear();
    arrival_time_distrb.reserve(number_photons);

    // distance < 25cm
    if (distance < 25) {
        // times are fixed shift i.e. direct path only
        double t_prop_correction = distance/vuv_vgroup_mean;
        for (int i = 0; i < number_photons; i++){
            arrival_time_distrb.push_back(t_prop_correction);
        }
    }
    // distance >= 25cm
    else {
        // determine nearest parameterisation in discretisation
        int index = std::round((distance - 25) / step_size);
        // check whether required parameterisation has been generated, generating if not
        if (VUV_timing[index].GetNdim() == 0) {
            generateparam(index);
        }
        // randomly sample parameterisation for each photon
        for (int i = 0; i < number_photons; i++){
            arrival_time_distrb.push_back(VUV_timing[index].GetRandom(VUV_min[index],VUV_max[index]));
        }  
    }
    return arrival_time_distrb;
}

// vis arrival times calculation function
vector<double> timeparam::getVisTime(TVector3 ScintPoint, TVector3 OpDetPoint, int number_photons) {
    
    // calculate point of reflection for shortest path accounting for difference in refractive indicies    
    // vectors for storing results
    TVector3 image(0,0,0);
    TVector3 bounce_point(0,0,0);
    TVector3 hotspot(0,0,0);
    TVector3 v_to_wall(0,0,0);
    
    // distance to wall    
    v_to_wall[0] = plane_depth - ScintPoint[0];

    // hotspot is point on wall where TPB is activated most intensely by the scintillation
    hotspot = ScintPoint + v_to_wall;
    
    // define "image" by reflecting over plane
    image = hotspot + v_to_wall*(n_LAr_vis/n_LAr_VUV);
    
    // find point of intersection with plane j of ray from the PMT to the image
    TVector3 tempvec = (OpDetPoint-image).Unit();
    double tempnorm= ((image-hotspot).Mag())/std::abs(tempvec[0]);
    bounce_point = image + tempvec*tempnorm;

    // calculate distance travelled by VUV light and by vis light
    double VUVdist = (bounce_point-ScintPoint).Mag();
    double Visdist = (OpDetPoint-bounce_point).Mag();

    // calculate times taken by each part
    vector<double> VUVTimes  = getVUVTime(VUVdist, number_photons);
    vector<double> ReflTimes(number_photons,0);
    double v;
    for (int i=0; i<number_photons; i++) {
        ReflTimes[i] = Visdist/vis_vmean;
    }

    // sum parts to get total transport times times                    
    vector<double> transport_time_vis(number_photons,0);
    for (int i=0; i<number_photons; i++) {
        transport_time_vis[i] = VUVTimes[i] + ReflTimes[i];
    }


    // *************************************************************************************************
    // Smearing of arrival time distribution
    // *************************************************************************************************
   
    gRandom->SetSeed(0);

    // calculate fastest time possible
    // vis part
    double vis_time = Visdist/vis_vmean;
    // vuv part
    double vuv_time;
    if (VUVdist < 25){
        vuv_time = VUVdist/vuv_vgroup_mean;
    }
    else {
        // find index of required parameterisation
        int index = std::round((VUVdist - 25) / step_size);
        // find shortest time
        vuv_time = VUV_min[index];
    }
    // sum
    double fastest_time = vis_time + vuv_time;

    // calculate angle between scintillation point and reflection point
    double delta_y = abs(ScintPoint[1] - bounce_point[1]);
    double delta_z = abs(ScintPoint[2] - bounce_point[2]);
    double delta_r = sqrt(pow(delta_y,2) + pow(delta_z,2));
    double theta = atan(delta_r/(plane_depth - ScintPoint[0])) * (180/3.1415);  //in degrees

    // calculate smearing parameters: 
    //  1). tau = exponential smearing factor
    //  2). delta = angular correction to range of random numbers generated, increases smearing as angle increases
    //  3). cutoff = largest smeared time allowed, preventing excessively large times caused by exponential
    double distance_cathode_plane = std::abs(plane_depth - ScintPoint[0]);
    double tau = 12.300 - 0.1144*distance_cathode_plane + 0.0004546*pow(distance_cathode_plane,2) - 6.1954e-7*pow(distance_cathode_plane,3);
    double delta = (0.003168 + 5.1967e-6*distance_cathode_plane + 2.205e-8*pow(distance_cathode_plane,2))*theta;
    double cutoff = 227.188 + 1.3925*distance_cathode_plane - 0.001521*pow(distance_cathode_plane,2);

    // apply smearing:
    for (int i=0; i < number_photons; i++){
        double arrival_time = transport_time_vis[i];
        double arrival_time_smeared;
        // if time is already greater than cutoff or minimum smeared time would be greater than cutoff, do not apply smearing
        if (arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(1.0-delta))-1) >= cutoff) {
            arrival_time_smeared = arrival_time;
        }
        // otherwise smear
        else {
            int counter = 0;
            // loop until time generated is within cutoff limit
            // most are within single attempt, very few take more than two, but this method could be made more efficient
            do {
                // don't attempt smearings too many times for cases near cutoff (very few cases)
                if (counter >= 10){
                    arrival_time_smeared = arrival_time; // don't smear
                    break;
                }
                else {
                    // generate random number in appropriate range            
                    double x = gRandom->Uniform(0.5,1.0 - delta);
                    // apply the exponential smearing
                    arrival_time_smeared = arrival_time + (arrival_time-fastest_time)*(exp(-tau*log(x))-1);
                }
                // increment counter
                counter++;
            }   while (arrival_time_smeared > cutoff);
        }
        transport_time_vis[i] = arrival_time_smeared;
    }  
    
    return transport_time_vis;
}