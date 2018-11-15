#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>

#include "utility_functions.h"
#include "library_access.h"
#include "libraryanalyze_light_histo.h"
#include "timeparam.h"


using namespace std;


// main driver for simulation


int main() { 
    
    //////////////////////////////////////////////////////////////////////////////
    ////////////------------ DEFINE REQUIRED FUNCTIONS ----------------///////////
    //////////////////////////////////////////////////////////////////////////////
    TF1 *fSpectrum = new TF1("fSpectrum",utility::SpectrumFunction,0,Q_Ar,1);//-----Beta decay spectrum
    TF1 *flandau_sn = new TF1("flandau_sn",utility::fsn, 0, 50, 1);//--SN Nu spectrum
    TF1 *flandau_so = new TF1("flandau_so",utility::fso, 0.2, 14.22, 1); // solar neutrino spectrum
    fSpectrum->SetParameter(0, Q_Ar);
    flandau_sn->SetParameter(0, Eav);
    flandau_so->SetParameter(0, Eav);

    TRandom3 *fGauss = new TRandom3();
    
    TF1 *fScintillation_function = new TF1("Scintillation Timing", utility::Scintillation_function, 0, scint_time_window, 3); 
    fScintillation_function->SetParameter(0, t_singlet); 
    fScintillation_function->SetParameter(1, t_triplet);    // t_singlet and t_triplet are defined in the header file (libraryanalyze_light_histo.h)
    if(gen_argon == true || fixed_energy == true || supernova == true || solar == true){ //i.e. an electron is the ionising particle
        fScintillation_function->FixParameter(2, 0); 
    }
    if(gen_radon == true){fScintillation_function->FixParameter(2, 1);} // i.e. an alpha particle is the ionising particle
    // the 3rd parameter is the type of particle 0 = electron, 1 = alpha (sets the ratio of the fast and slow components - see utility_functions.cc)


    //////////////////////////////////////////////////////////////////////////////
    ////////////------INTRODUCE SIMULATION/GENERATION CONDITIONS-------///////////
    //////////////////////////////////////////////////////////////////////////////
    int max_events;
    int scint_yield;
    string particle;

    // OUTPUT what is being simulationed
    if(fixed_energy == true) {
        max_events = max_events_FE;
        scint_yield = scint_yield_electron;
        particle = "electron";
        cout << endl << "Generating " << max_events << " events, of fixed energy: " << fixedE << " MeV." << endl;
    }
    if(gen_argon == true) {
        max_events = max_events_Ar;
        scint_yield = scint_yield_electron;
        particle = "electron";
        cout << endl << "Generating " << max_events << ", Ar 39 decays in time window: " << time_window << " seconds." << endl;
        cout << "This is equal to " << time_frames << " PMT readout frames." << endl;
        cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
        cout << "***NOTE. In ONE TPC, We expect to see: " << Ar_decays_per_sec << " Ar 39 decays each second.***" << endl;
        cout << "//////////////////////////////////////////////////////////////////////////////" << endl << endl;
    }
    if(supernova == true){
        max_events = max_events_SN;
        scint_yield = scint_yield_electron;
        particle = "electron";
        cout << "\nGenerating " << max_events << ", supernova events.\n";
    }
    if(solar == true){
        max_events = max_events_SO;
        scint_yield = scint_yield_electron;
        particle = "electron";
        cout << "\nGenerating " << max_events << ", solar events.\n";
    }
    if(gen_radon == true){ //gen_radon == true
        max_events = max_events_Rn;
        scint_yield = scint_yield_alpha;
        particle = "alpha";
        cout << "\nGenerating " << max_events << ", Radon 222 decays in time window: " << time_window << " seconds." << endl;
        cout << "This is equal to " << time_frames << " PMT readout frames." << endl;
        cout << endl << "/////////////////////////////////////////////////////////////////////////////" << endl;
        cout << "***NOTE. In ONE TPC, we expect to see: " << Rn_decays_per_sec << " radon-222 decays each second.***" << endl;
        cout << "/////////////////////////////////////////////////////////////////////////////" << endl << endl;
    }

    // OUTPUT the foil configuration chosen
    cout << "With the following set up: ";
    if(config == 0) {cout << "Full Foils" << endl; }
    if(config == 1) {cout << "Cath Foils" << endl; }
    if(config == 2) {cout << "VUV only" << endl; }
    cout << endl;
    if (config == 0 || config ==1) {
        reflected = true;
        reflT = true;
    }
    if (config ==2) {
        reflected = false;
        reflT = false;
    }


    ////////////////////////////////////////////////////////////////////////////
    ////////////--------------CREATING THE TREES-----------------///////////////
    ////////////////////////////////////////////////////////////////////////////
    // IMPORTANT \\
    // event_tree will store an entry for each event (i.e. a supernova neutrino)
    // data_tree will store an entry for each detected photon from ALL events simulated 
    // => data_tree is MUCH larger as each event will produce ~24,000 photons/MeV (but only a fraction of these are "detected")
    // It is the optical libraries which determine how many photons are detected by each PMT.

    // data_tree was designed such that each photoelectron would have their own entry in the Ttree.
    // each photoelectron entry includes: the time it took (scint time + transport time), which pmt it hit, (x,y,z) origin, and which event it came from (data_event)    A current issue with the event number tracking is that it does not carry over between several iterations of the code. For single uses this is fine, but memory consuption and speed for running this code a single time becomes a limiting factor, such that running many instances of this code    //in parallel becomes much more useful. Therefore carrying over the event    //number between iterations may be useful.

    // ------ DATA TREES ------
    ////////////////////////////
    // data tree for TOTAL (VUV + Vis) light from ALL events generated
    ///////////////////////////
    if (do_timings) { data_tree->Branch("data_time", &data_time, "data_time/D"); }
    data_tree->Branch("data_pmt", &data_pmt, "data_pmt/I");
    data_tree->Branch("data_x_pos", &data_x_pos, "data_x_pos/D");
    data_tree->Branch("data_y_pos", &data_y_pos, "data_y_pos/D");
    data_tree->Branch("data_z_pos", &data_z_pos, "data_z_pos/D"); 
    data_tree->Branch("data_event", &data_event, "data_event/I");

    /////////////////////////
    // datatree for VUV light from ALL events generated
    ////////////////////////// 
    if (do_timings) { data_tree_vuv->Branch("data_time_vuv", &data_time_vuv, "data_time_vuv/D"); }
    data_tree_vuv->Branch("data_pmt_vuv", &data_pmt_vuv, "data_pmt_vuv/I");
    data_tree_vuv->Branch("data_x_pos_vuv", &data_x_pos_vuv, "data_x_pos_vuv/D"); 
    data_tree_vuv->Branch("data_y_pos_vuv", &data_y_pos_vuv, "data_y_pos_vuv/D");
    data_tree_vuv->Branch("data_z_pos_vuv", &data_z_pos_vuv, "data_z_pos_vuv/D");
    data_tree_vuv->Branch("data_event_vuv", &data_event_vuv, "data_event_vuv/I");

    //////////////////////////////////
    // datatree for visible light from ALL events generated
    /////////////////////////////////
    if (do_timings) { data_tree_vis->Branch("data_time_vis", &data_time_vis, "data_time_vis/D"); }
    data_tree_vis->Branch("data_pmt_vis", &data_pmt_vis, "data_pmt_vis/I");
    data_tree_vis->Branch("data_x_pos_vis", &data_x_pos_vis, "data_x_pos_vis/D");
    data_tree_vis->Branch("data_y_pos_vis", &data_y_pos_vis, "data_y_pos_vis/D");
    data_tree_vis->Branch("data_z_pos_vis", &data_z_pos_vis, "data_z_pos_vis/D");
    data_tree_vis->Branch("data_event_vis", &data_event_vis, "data_event_vis/I");

    // ------------ EVENT TREE -------------
    //event_tree is designed to save one entry for every event (e.g. a supernova event).
    // Each event stores: the event number, the voxel number which the event occured in, the (x,y,z) positions of the event and the energy of the event.
    event_tree->Branch("event_no", &event_no, "event_no/I");
    event_tree->Branch("event_vox", &event_vox, "event_vox/I");
    event_tree->Branch("event_x_pos", &event_x_pos, "event_x_pos/D");
    event_tree->Branch("event_y_pos", &event_y_pos, "event_y_pos/D");
    event_tree->Branch("event_z_pos", &event_z_pos, "event_z_pos/D");
    event_tree->Branch("event_E", &event_E, "event_E/D");


    ////////////////////////////////////////////////////////////////////////////////////////
    ////////////------------lOADING THE DESIRED OPTICAL LIBRARY------------------///////////
    ////////////////////////////////////////////////////////////////////////////////////////
    if(config == 0) {libraryfile = "dune_withfoils_lib_cleaned.root"; }
    if(config == 1) {libraryfile = "dune_withfoils_lib_cleaned.root"; }          
    if(config == 2) {libraryfile = "dune_withfoils_lib_cleaned.root"; }
    lar_light.LoadLibraryFromFile(libraryfile, reflected, reflT);


    if (do_timings) {
        ////////////////////////////////////////////////////////////////////////////////////
        ////////////-------------FILLING A VECTOR OF PMT POSITIONS---------------///////////
        ////////////////////////////////////////////////////////////////////////////////////
        //NOTE: the .txt file is just a list of the PMT positions, which will simply be filled into a vector.
        cout << "Loading Photon Detector positions...\n";
        ifstream myfile;
        myfile.open("bars_positions_dune1x2x6.txt"); // ALTERED FOR DUNE
        if(myfile.is_open()) {cout << "File opened successfully" << endl; }
        while(!myfile.eof()) {
            double num_pmt, x_pmt, y_pmt, z_pmt;
            if(myfile >> num_pmt >> x_pmt >> y_pmt >> z_pmt) {
                vector<double> line_data({num_pmt, x_pmt, y_pmt, z_pmt});
                myfile_data.push_back(line_data);
            }
                else{break; }
        }
        myfile.close();
        cout << "Positions Loaded!\n\n";
        // End Reading out positions of PMT from txt file.
    }

    // initalise timings class
    timeparam Timing(step_size);


    //////////////////////////////////////////////////////////////////////////////
    ////////////-------------------MAIN CODE---------------------------///////////
    //////////////////////////////////////////////////////////////////////////////
    cout << "Generating event positions and energies:\n";
    gRandom->SetSeed(0);


    ///////// -----------------  Create events --------------------- ///////////

    energy_list.reserve(max_events);
    if (do_timings) { decay_time_list.reserve(max_events); }
    voxel_list.reserve(max_events);
    position_list.resize(max_events,vector<double>(3,0.));
    // these vectors were defined in the header file (libraryanalyze_light_histo.h)

    // a loop to deal with each individual event: basically gets the details of the event and puts them into the event_tree
    int j = 1;
    for(int event = 0; event < max_events; event++) {
        
        // output completion %
        if ( (event != 0) && (max_events >= 10) &&  (event % (max_events/10) == 0) ) {
            cout << Form("%i0%% Completed...\n", j);
            j++;
        }
        
        // DETERMINE THE ENERGY OF THE EVENT
        double energy;
        if(fixed_energy == true) {energy = fixedE;} 
        if(gen_argon == true) {energy = fSpectrum->GetRandom();}    // pull from the Ar beta spectrum (see utility_functions.cc)
        if(supernova == true) {energy = flandau_sn->GetRandom();}   // pull from the predicted SN spectrum (see utility_functions.cc)
        if(solar == true) {energy = flandau_so->GetRandom();}		// pull from solar neutrino spectrum
        if(gen_radon == true) {energy = fGauss->Gaus(Q_Rn, 0.05);}  // Gaus(av,sigma) - is a ROOT function, pulls from a Gaussian
         

        // DETERMINE THE POSITION OF THE VOXEL IN WHICH THE EVENT OCCURRED & THE VOXEL NUMBER 
        vector<double> position(3); //creates an empty vector of length 3
        int rand_voxel;

        // 3 possible cases: random (x,y,z), fixed x & random (y,z) and fixed (x,y,z) - this choice is made in the header file
        if(random_pos == true) {
            // choose random position in 3D range
            for (int i=0; i<3; i++) {
                position[i] = gRandom->Uniform(PosMin[i],PosMax[i]);
            }
        }
        else if(fixed_xpos == true) { 
            // choose a random voxel with a fixed x (drift distance) position.
            double randomY = ((rand() % 400) - 200)+0.5;    // random Y voxel
            double randomZ = (rand() % 500) + 0.5;          // random Z voxel
            position[0] = fixedX; position[1]= randomY; position[2] = randomZ; // fill the array
        }
        else { // fixed_pos == true
        position[0] = fixedX; position[1]= fixedY; position[2] = fixedZ;
        }
        
        double pos_ar[3] = {position[0],position[1],position[2]};
        rand_voxel = lar_light.GetVoxelID(pos_ar); // get the ID of the voxel

        // fill the vectors
        position_list[event] = position;
        voxel_list.push_back(rand_voxel); // push back the position of the voxel onto an array
        energy_list.push_back(energy); // push back the energy of the event onto the array


        if (do_timings) {
        //double decay_time = time_window * gRandom->Uniform(1.); // selects a random time of the event to decay within the specified time window
             double decay_time = 0; // all decays happen at t = 0 (obviously this is not realistic)
             decay_time_list.push_back(decay_time);
        }

        //////////////////////////////
        ///// Fill the event tree ////
        //////////////////////////////
        event_no = event; // event number is equal to the current loop iteration
        event_vox = rand_voxel;
        event_x_pos = position[0];
        event_y_pos = position[1];
        event_z_pos = position[2];
        event_E = energy_list.at(event); // pull from the energy list at the current iteration
        event_tree->Fill();

    } // end of event loop 


    ////////*********************************************************************************************************************************
    ///////////////////////////////////////////
    ////////NOW THE DATA TREE WILL BE FILLED///
    ///////////////////////////////////////////

    // What is going to happen below:
    // We are going to take an event, and then go the first PMT and see how many visible & VUV photons from THAT event hit THAT PMT
    // We will then go onto the next PMT and do the same
    // and then the next PMT... and so on.
    // We will then do all of this for the next event, and then the next event, and so on...
    
    //FOR EFFICIENCY DEBUGGING: SETUP CLOCKS  
    chrono::steady_clock::time_point t_all_i;
    chrono::steady_clock::time_point t_all_f;
    chrono::steady_clock::time_point t_VUVTime_i;
    chrono::steady_clock::time_point t_VUVTime_f;
    chrono::steady_clock::time_point t_ReflTime_i;
    chrono::steady_clock::time_point t_ReflTime_f;
    chrono::duration<double> timespan_all (0.);
    chrono::duration<double> timespan_VUV (0.);
    chrono::duration<double> timespan_Refl (0.);
    
    t_all_i = chrono::steady_clock::now();
    
    if (do_timings) { 
        cout << "Determining number of photon hits and timings:\n";
    }
    else            { cout << "Determining number of photon hits:\n"; }

    //Loop over each PMT for each event
    j = 1;    
    for(int events = 0; events < max_events; events++) {
        
        if ( (events != 0) && (max_events >= 10) &&  (events % (max_events/10) == 0) ) {
            cout << Form("%i0%% Completed...\n", j);
            j++;
        }
        
        // Begin looping over the PMT array.
        for(int pmt_loop = 0; pmt_loop < 120; pmt_loop++) { //  120 optical channels for dune bars configuration

            int num_pmt = pmt_loop;

            // - This function (defined in library_access.cc) will determine how many VUV and Visble photons hit the given PMT
            vector<double> pmt_hits = lar_light.PhotonLibraryAnalyzer(energy_list.at(events), scint_yield, quantum_efficiency, catcov, vuvfrac, visfrac, num_pmt, voxel_list.at(events));
            // NB this vector has form: (no of VUV photons, no of visible photons)

            int num_VUV = pmt_hits.at(0);
            int num_VIS = pmt_hits.at(1);

            // if no photons from this event for this PMT, go to the next event.
            if(num_VUV+num_VIS == 0) {continue; } // forces the next iteration
            
            // position and event number
            data_x_pos = position_list[events][0];
            data_x_pos_vuv = position_list[events][0];
            data_x_pos_vis = position_list[events][0];

            data_y_pos = position_list[events][1];
            data_y_pos_vuv = position_list[events][1];
            data_y_pos_vis = position_list[events][1];

            data_z_pos = position_list[events][2];
            data_z_pos_vuv = position_list[events][2];
            data_z_pos_vis = position_list[events][2];

            data_event = events;
            data_event_vuv = events;
            data_event_vis = events;


            if (do_timings) {
                ///*************************** ///
                ///***** TIMING OPERATIONS *** ///
                ///*************************** ///
                // NOTE:    Distances in cm and times in ns, returns a vector of doubles for VUV arrival times OF EACH PHOTON
                // Get the (x,y,z) position of the scintillation point and of PMT as we need this to work out transport time
                TVector3 ScintPoint(data_x_pos,data_y_pos,data_z_pos);
                TVector3 OpDetPoint(myfile_data.at(num_pmt).at(1),myfile_data.at(num_pmt).at(2),myfile_data.at(num_pmt).at(3));
                
                ///////////////////////
                // VUV photons times
                ///////////////////////
                
                // transport times
                double distance_to_pmt = (OpDetPoint-ScintPoint).Mag();   
                vector<double> transport_time_vuv;
                if(num_VUV != 0) {
                    t_VUVTime_i = chrono::steady_clock::now();
                    transport_time_vuv = Timing.getVUVTime(distance_to_pmt, num_VUV);
                    t_VUVTime_f = chrono::steady_clock::now();
                    timespan_VUV += chrono::duration_cast<chrono::duration<double>>(t_VUVTime_f-t_VUVTime_i);
                }
                
                // total times
                double total_time_vuv; // total time for each VUV photon to hit a pmt from creation
                for(auto& x: transport_time_vuv) { //we're just looping through the transport_time_vuv vector here

                    // data_time has several parts and is simply: transport time + decay time + scintillation time
                    // x is the transport time (in NANOSECONDS) - x * 0.001 converts from ns -> micros_s
                    // the scintillation function timing is also converted to microseconds
                    // the time window offset is when the decay occured, given already in microseconds            
                    total_time_vuv = (x*0.001+(decay_time_list.at(events) + fScintillation_function->GetRandom())*1000000.); // in microseconds 

                    // apply time cut if required
                    // 0.1 microseconds = 100 ns! 
                    if(cut == true && total_time_vuv > time_cut){ continue; } // continue to the next iteration without filling - don't bother filling it
                    
                    // fill data trees
                    data_time = total_time_vuv;
                    data_time_vuv = total_time_vuv;
                    data_pmt = num_pmt;
                    data_pmt_vuv = num_pmt;
                    data_tree_vuv->Fill();
                    data_tree->Fill();
                } //end of looping over the transport time vector
                transport_time_vuv.clear();


                ////////////////////////
                // visible photons times 
                ////////////////////////

                vector<double> transport_time_vis;
                if(num_VIS != 0 && config == 1) {  // NOTE config == 1 is cathode foils configuration
                    
                    // calculate the transport time
                    t_ReflTime_i = chrono::steady_clock::now();
                    vector<double> transport_time_vis = Timing.getVisTime(ScintPoint, OpDetPoint, num_VIS);
                    t_ReflTime_f = chrono::steady_clock::now();
                    timespan_Refl += chrono::duration_cast<chrono::duration<double>>(t_ReflTime_f-t_ReflTime_i);
                    
                    // calculate the total time
                    double total_time_vis;
                    for(auto &y : transport_time_vis) { //looping through the transport_time_vis vector

                    	// as with VUV times, calculate total time
                        total_time_vis = (y*0.001+(decay_time_list.at(events) + fScintillation_function->GetRandom())*1000000.); // in microseconds

                        // apply time cut if required 
                        if(cut == true && total_time_vis > time_cut){ continue; } // go onto the next interation - cut has been made
                        
                        // fill data trees
                        data_time = total_time_vis; 
                        data_time_vis = total_time_vis;
                        data_pmt = num_pmt;
                        data_pmt_vis = num_pmt;
                        data_tree_vis->Fill();
                        data_tree->Fill();
                    } // end of loop through transport_time_vis vector
                    transport_time_vis.clear();
                }
                
            }
            else { // not doing timing calculations
                // fill data trees
                data_pmt = num_pmt;
                data_pmt_vuv = num_pmt;
                data_pmt_vis = num_pmt;
                
                // adds position and event data to appropriate data_trees for each photon detected  
		        if(num_VUV != 0) {
			         for (int i = 0; i < num_VUV; i++){
                        data_tree_vuv->Fill();
                        data_tree->Fill();   
                     }
		        }
		        if(num_VIS != 0) {
                	 for (int i = 0; i < num_VIS; i++){
                        data_tree_vis->Fill();
                        data_tree->Fill();   
                     }
		        }
                
            }
        
        } //end pmt loop
    } //end loop over events
    

    // display time taken - for debugging purposes
    t_all_f = std::chrono::steady_clock::now();
    timespan_all = chrono::duration_cast<chrono::duration<double>>(t_all_f-t_all_i);
    cout << "Total data tree filling took " << timespan_all.count() << " seconds\n";
    cout << "Direct light timings calculation took " << timespan_VUV.count() << " seconds in total\n";
    cout << "Reflected light timings calculation took " << timespan_Refl.count() << " seconds in total\n";
    
    energy_list.clear();
    voxel_list.clear();
    if (do_timings) { decay_time_list.clear(); }
    
    // write the files
    event_file.Write();

    return 0;

} // end main
