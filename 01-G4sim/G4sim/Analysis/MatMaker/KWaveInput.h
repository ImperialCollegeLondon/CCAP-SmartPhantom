#ifndef KWAVEINPUT_H
#define KWAVEINPUT_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <ctime>

#include "hdf5_hl.h"
#include "fftw3.h"
#include "Linterp/linterp.h"

// Class to handle writing k-Wave input
class KWaveInput
{
public:
    KWaveInput();
    // Create + open .h5 file
    explicit KWaveInput(const std::string& fn);
    ~KWaveInput();
    
    void WriteInput();

    
    void circshift(double *out, const double *in, int xdim, int ydim, int xshift, int yshift);
    void circshift(double *out, const double *in, int xdim, int ydim, int zdim, int xshift, int yshift, int zshift);

    void fftshift(double *out, const double *in, int xdim, int ydim);
    void fftshift(double *out, const double *in, int xdim, int ydim, int zdim);

    void ifftshift(double *out, const double *in, int xdim, int ydim);
    void ifftshift(double *out, const double *in, int xdim, int ydim, int zdim);
    
    std::vector<double> Linspace(double first, double last, int len);

    double* Blackman(int wSize);
    void WindowFilter(int* nArr, fftw_complex* data);
    void RotWindowFilter(int* nArr, fftw_complex* data);
    float* Smooth(double* energySum, int* binDim, double grueneisen = 1, bool restoreMagnitude = true, bool rotWindow = true);
    
    void WriteHeader();
    
    // Calculating variables
    // Calculate dt
    float CalculateDT(float cfl, float* gridSpacing, float speed);

    // Calculate Nt
    unsigned long CalculateNT(float tend, float dt);
    
    // For calculating sensor indices (assuming matlab indexing for input which starts from 1)
    // A rough implementation of sensor.mask() and find(sensor.mask) in matlab
    std::vector<unsigned long> SensorMaskIndices(unsigned long nx, int senMinX, int senMaxX, int senDx,
                                                 unsigned long ny, int senMinY, int senMaxY, int senDy,
                                                 unsigned long nz, int senMinZ, int senMaxZ, int senDz);
    
    // For writing k-wave 1D float variable
    void Write1DFloat(std::string datasetName, float val);

    // For writing k-wave 1D long variable
    void Write1DLong(std::string datasetName, unsigned long val);

    // For writing k-wave mult dim long variable (i.e. sensor_mask_index)
    void WriteMultDLong(std::string datasetName, std::vector<unsigned long>& indicesVec);
    
    // For writing data (i.e. p0_source_input)
    void WritePData(std::string datasetName, int nx, int ny, int nz, float* data);
    
    // For writing file
    void WriteFile();
    
    void SetBinDim(int* dim) { binDim = dim; };
    void SetEnergyData(double* data) { energyData = data; };
    void SetSmooth(bool val) { enableSmooth = val; };
    void SetGrueneisen(double val) { grueneisen = val; };
    void SetRestoreMagnitude(bool val) { restoreMagnitude = val; };
    void SetRotWindow(bool val) { rotWindow = val; };
    void SetNx(unsigned long val) { Nx = val; };
    void SetNy(unsigned long val) { Ny = val; };
    void SetNz(unsigned long val) { Nz = val; };
    void SetDx(float val) { dx = val; };
    void SetDy(float val) { dy = val; };
    void SetDz(float val) { dz = val; };
    void SetPMLSizeX(unsigned long val) { pmlSizeX = val; };
    void SetPMLSizeY(unsigned long val) { pmlSizeY = val; };
    void SetPMLSizeZ(unsigned long val) { pmlSizeZ = val; };
    void SetPMLAlphaX(float val) { pmlAlphaX = val; };
    void SetPMLAlphaY(float val) { pmlAlphaY = val; };
    void SetPMLAlphaZ(float val) { pmlAlphaZ = val; };
    void SetMediumSoundSpeed(float val) { mediumSoundSpeed = val; };
    void SetCRef(float val) { cRef = val; };
    void SetTEnd(float val) { tEnd = val; };
    void SetCFL(float val) { cfl = val; };
    void SetSensorMask(std::vector<unsigned long>& val) { sensorMask = val; };
        
    unsigned long GetNx() { return Nx; };
    unsigned long GetNy() { return Ny; };
    unsigned long GetNz() { return Nz; };
    float GetDx() { return dx; };
    float GetDy() { return dy; };
    float GetDz() { return dz; };
    unsigned long GetPMLSizeX() { return pmlSizeX; };
    unsigned long GetPMLSizeY() { return pmlSizeY; };
    unsigned long GetPMLSizeZ() { return pmlSizeZ; };
    float GetPMLAlphaX() { return pmlAlphaX; };
    float GetPMLAlphaY() { return pmlAlphaY; };
    float GetPMLAlphaZ() { return pmlAlphaZ; };
    float GetMediumSoundSpeed() { return mediumSoundSpeed; };
    float GetCRef() { return cRef; };
    float GetTEnd() { return tEnd; };
    float GetCFL() { return cfl; };
    float GetDt() { return dt; };
    unsigned long GetNt() { return Nt; };
    
private:
    hid_t fileId;
    
    int* binDim;
    double* energyData;
    float* fenergyData;
    
    bool enableSmooth;
    double grueneisen;
    bool restoreMagnitude;
    bool rotWindow;
    
    unsigned long Nx;
    unsigned long Ny;
    unsigned long Nz;
    
    float dx;
    float dy;
    float dz;
    
    unsigned long pmlSizeX;
    unsigned long pmlSizeY;
    unsigned long pmlSizeZ;
    float pmlAlphaX;
    float pmlAlphaY;
    float pmlAlphaZ;
    
    float mediumSoundSpeed;
    float cRef;
    
    float tEnd;
    float cfl;
    
    float dt;
    unsigned long Nt;
    
    std::vector<unsigned long> sensorMask;
};
#endif // KWAVEINPUT_H
