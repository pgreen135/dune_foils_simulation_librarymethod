#ifndef LIBRARY_ACCESS_H
#define LIBRARY_ACCESS_H

#include <string>
#include <vector>

//This file is designed to access the visibility parameters from the
//optical libraries, which are needed to calculate the number of photoelectrons
//incident on each PMT.

class LibraryAccess{

  public:
    void LoadLibraryFromFile(std::string libraryfile, bool reflected, bool reflT0);
    const float* GetReflT0(size_t Voxel, int no_pmt);
    const float* GetReflCounts(size_t Voxel, int no_pmt, bool is_reflT0);
    const float* GetCounts(size_t Voxel, int no_pmt);
    const float* GetLibraryEntries(int VoxID, bool wantReflected, int no_pmt);
    std::vector<int> GetVoxelCoords(int id, double position[3]);
    int GetVoxelID(double* Position);
    std::vector<double> PhotonLibraryAnalyzer(double _energy, const int _scint_yield, const double _quantum_efficiency, const double _catcov, const double _vuvfrac, const double _visfrac, int _pmt_number, int _rand_voxel);

    LibraryAccess();

  private:
    std::vector<std::vector<float> > table_;
    std::vector<std::vector<float> > reflected_table_;
    std::vector<std::vector<float> > reflT_table_;
    
    // dune optical library geometry
    const double gLowerCorner[3] = {0., -658.099, -302.946};
    const double gUpperCorner[3] = {379.662, 699.599, 1443.53};
    const int gxSteps = 75;
    const int gySteps = 125;
    const int gzSteps = 300;

};

#endif
