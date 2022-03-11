#include "KWaveInput.h"

KWaveInput::KWaveInput()
{}

KWaveInput::KWaveInput(const std::string& fn)
{
    fileId = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    WriteHeader();
}

KWaveInput::~KWaveInput()
{
    std::cout << "Beginning destructor for KWaveInput." << std::endl;
    delete fenergyData;
    H5Fclose(fileId);
    std::cout << "End of destructor for KWaveInput." << std::endl;    
}


// Taken from: 
// https://stackoverflow.com/questions/5915125/fftshift-ifftshift-c-c-source-code

void KWaveInput::circshift(double *out, const double *in, int xdim, int ydim, int xshift, int yshift)
{
    for (int i = 0; i < xdim; i++) 
    {
        int ii = (i + xshift) % xdim;
        for (int j = 0; j < ydim; j++) 
        {
        int jj = (j + yshift) % ydim;
        out[ii * ydim + jj] = in[i * ydim + j];
        }
    }
}

void KWaveInput::circshift(double *out, const double *in, int xdim, int ydim, int zdim, int xshift, int yshift, int zshift)
{
    for (int i = 0; i < xdim; i++) 
    {
        int ii = (i + xshift) % xdim;
        for (int j = 0; j < ydim; j++) 
        {
            int jj = (j + yshift) % ydim;
            for(int k=0; k< zdim; k++)
            {
                int kk = (k + zshift) % zdim;
                int index2 = (kk + zdim*(jj+ydim*ii));
                int index = k + zdim*(j+ydim*i);
                out[index2] = in[index];
            }
        }
    }
}

void KWaveInput::fftshift(double *out, const double *in, int xdim, int ydim)
{
    return circshift(out, in, xdim, ydim, (xdim/2), (ydim/2));
}

void KWaveInput::fftshift(double *out, const double *in, int xdim, int ydim, int zdim)
{
    return circshift(out, in, xdim, ydim, zdim, (xdim/2), (ydim/2), (zdim/2));
}

void KWaveInput::ifftshift(double *out, const double *in, int xdim, int ydim)
{
    return circshift(out, in, xdim, ydim, ((xdim+1)/2), ((ydim+1)/2));
}

void KWaveInput::ifftshift(double *out, const double *in, int xdim, int ydim, int zdim)
{
    return circshift(out, in, xdim, ydim, zdim, ((xdim+1)/2), ((ydim+1)/2), ((zdim+1)/2));
}

// Return an evenly spaced 1-d grid of doubles
// taken from https://rncarpio.github.io/linterp/
std::vector<double> KWaveInput::Linspace(double first, double last, int len)
{
    std::vector<double> result(len);
    double step = (last - first) / (len - 1);
    for(int i=0; i<len; i++)
    {
        result[i] = first + i*step;
    }
    return result;
}

// Meant to produce the same result as blackman from scipy.signal, producing symmetric window
double* KWaveInput::Blackman(int wSize)
{
    int windowSize = wSize; // Number of points in output window
    double* window = new double[windowSize];
    
    if( windowSize%2 == 0) // If even, return periodic window
    {
        windowSize += 1;
    }
    
    double stepsize = (double)windowSize/(windowSize - 1);
    const double pi = std::acos(-1);
    
    double a0 = 0.42;
    double a1 = 0.5;
    double a2 = 0.08;
    
    double n = 0;
    double maxVal = 0;
    
    for(int i=0; i<wSize; i++) // counter is wrt wSize to handle even and odd as specified in smooth.m from k-wave
    {
        double blackman = a0 - a1*std::cos(2*pi*n/windowSize) + a2*std::cos(4*pi*n/windowSize);
        //blackman += std::numeric_limits<float>::epsilon(); // machine epsilon for float to make values positive
                                                           // copying python example of: kwave_data_filters.py
        window[i] = blackman;
        n += stepsize;
        if(maxVal < blackman)
            maxVal = blackman;
    }

    return window;
}

void KWaveInput::WindowFilter(int* nArr, fftw_complex* data)
{
    // To create Non-symmetrical window for even dimensions
    int dim = sizeof(nArr)/sizeof(nArr[0]) + 1; // calculate number of dimensions, 
                                                // though for now code only works for 3d

    int adjustData[dim];
    for(int i=0; i<dim; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }

    double filterX[adjustData[0]] = {0};
    double filterY[adjustData[1]] = {0};
    double filterZ[adjustData[2]] = {0};

    for(int axis=0; axis<dim; axis++)
    {
        if(axis==0)
        {
            double* blackman = Blackman(adjustData[0]);
            //double ifftBlackman[windowSize] = {0};
            //ifftshift(ifftBlackman, blackman, windowSize, 1);
            
            for(int i=0; i<adjustData[0]; i++)
            {
                filterX[i] = blackman[i];
            }
        }
        else if(axis==1)
        {
            double* blackman = Blackman(adjustData[1]);
            //double ifftBlackman[windowSize] = {0};
            //ifftshift(ifftBlackman, blackman, windowSize, 1);
            
            for(int i=0; i<adjustData[1]; i++)
            {
                filterY[i] = blackman[i];
            }
        }
        else if(axis==2)
        {
            double* blackman = Blackman(adjustData[2]);
            //double ifftBlackman[windowSize] = {0};
            //ifftshift(ifftBlackman, blackman, windowSize, 1);
            
            for(int i=0; i<adjustData[2]; i++)
            {
                filterZ[i] = blackman[i];
            }
        }
        else
        {
            std::cout << "I can't handle more or less than three dimensions yet!" << std::endl;
        }
    }
    
    double filterCombined[adjustData[0]*adjustData[1]*adjustData[2]];
    for(int i=0; i<adjustData[0]; i++)
    {
        for(int j=0; j<adjustData[1]; j++)
        {
            for(int k=0; k<adjustData[2]; k++)
            {
                int index = k + adjustData[2]*(j+adjustData[1]*i);
                filterCombined[index] = std::abs(filterX[i]*filterY[j]*filterZ[k]);
                //std::cout << i << '\t' << j << '\t' << k << '\t' << filterCombined[index] << std::endl;
            }
        }
    }
    
    // Messy way of slicing filterCombined to input dimensions
    double filterSliced[nArr[0]*nArr[1]*nArr[2]];
    int sliceIndex = 0;

    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int index = k + nArr[2]*(j+nArr[1]*i);
                filterSliced[index] = filterCombined[sliceIndex];
                //std::cout << i << '\t' << j << '\t' << k << '\t' << filterSliced[index] << std::endl;
                sliceIndex++;
            }
            //std::cout << std::endl;
            if(nArr[2]!=adjustData[2])
                sliceIndex++;
        }
        //std::cout << std::endl;
        if(nArr[1]!=adjustData[1])
            sliceIndex += adjustData[2];
    }

    double ifftFilterCombined[nArr[0]*nArr[1]*nArr[2]];
    ifftshift(ifftFilterCombined, filterSliced, nArr[0], nArr[1], nArr[2]);    
    
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int index = k+nArr[2]*(j+nArr[1]*i);
                data[index][0] *= ifftFilterCombined[index]; // Real
                data[index][1] *= ifftFilterCombined[index]; // Imaginary
                //std::cout << i << '\t' << j << '\t' << k << '\t' << data[index][0] << '\t' << data[index][1] << std::endl;
            }
        }
    }
}

void KWaveInput::RotWindowFilter(int* nArr, fftw_complex* data)
{    
    // To create Non-symmetrical window for even dimensions
    int dim = sizeof(nArr)/sizeof(nArr[0]) + 1; // calculate number of dimensions, 
                                                // though for now code only works for 3d

    // Adjust dimensions to account for nonsymmetric filter if dimension has even number of terms
    int adjustData[dim];
    //std::cout << nArr[0] << '\t' << nArr[1] << '\t' << nArr[2] << std::endl;
    for(int i=0; i<dim; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }

    int L = *std::max_element(adjustData,adjustData+dim);
    double radius = (L-1)/2.;
    std::vector<double> ll = Linspace(-radius,radius,L);
    std::vector<double> xx = Linspace(-radius,radius,adjustData[0]);
    std::vector<double> yy = Linspace(-radius,radius,adjustData[1]);
    std::vector<double> zz = Linspace(-radius,radius,adjustData[2]);

    std::vector<double> r;
    r.reserve(adjustData[0]*adjustData[1]*adjustData[2]);

    // Create Blackman window and put into vector winLin 
    std::vector<double> winLin;
    winLin.reserve(L);
    double* bWin = Blackman(L);

    for(int i=0; i<L; i++)
    {
        winLin[i] = bWin[i];
        //std::cout << winLin[i] << std::endl;
    }
    
    // Set up grid based on radius
    for(int i=0; i<adjustData[0]; i++)
    {
        for(int j=0; j<adjustData[1]; j++)
        {
            for(int k=0; k<adjustData[2]; k++)
            {
                int index = k+adjustData[2]*(j+adjustData[1]*i);
                double rr = std::sqrt( xx[i]*xx[i] + yy[j]*yy[j] + zz[k]*zz[k] );
                if( rr > radius)
                    rr = radius;
                r[index] = rr;
                //std::cout << rr << '\t';
            }
            //std::cout << std::endl;
        }
        //std::cout << std::endl;
    }
    
    // Set up interpolator parameters
    std::vector<std::vector<double>::iterator> gridIterList;
    gridIterList.push_back(ll.begin());
    array<int,1> gridSizes;
    gridSizes[0] = L;
    InterpMultilinear<1, double> 
    interp_ML(gridIterList.begin(),gridSizes.begin(),winLin.data(),winLin.data()+L);

    // Messy way of slicing 'r' dimensions to 'nArr'
    double rotWin[nArr[0]*nArr[1]*nArr[2]];
    int rIndex = 0;
    //int rotIndex = 0;

    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int rotIndex = k + nArr[2]*(j+nArr[1]*i);
                array<double,1> args = {r[rIndex]};
                //std::cout << r[rIndex] << '\t';
                rotWin[rotIndex] = std::abs(interp_ML.interp(args.begin()));
                //std::cout << i << '\t' << j << '\t' << k << '\t' << rotWin[rotIndex] << std::endl;
                //std::cout << rotIndex << '\t' << rotWin[rotIndex] << '\t';
                //std::cout << rotIndex << '\t';
                rIndex++;
            }
            //std::cout << std::endl;
            if(nArr[2]!=adjustData[2])
                rIndex++;
        }
        //std::cout << std::endl;
        if(nArr[1]!=adjustData[1])
            rIndex += adjustData[2];
    }

    // Perform ifftshift of window filter
    double rotWin2[nArr[0]*nArr[1]*nArr[2]] = {0};
    ifftshift(rotWin2, rotWin, nArr[0], nArr[1], nArr[2]);    
    
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int index = k+nArr[2]*(j+nArr[1]*i);
                //std::cout << index << '\t';
                //std::cout << rotWin2[index] << '\t';
                data[index][0] *= rotWin2[index]; // Real
                data[index][1] *= rotWin2[index]; // Imaginary
            }
            //std::cout << std::endl;
        }
        //std::cout << std::endl;
    }
}

float* KWaveInput::Smooth(double* energySum, int* binDim, double grueneisen, bool restoreMagnitude, bool rotWindow)
{
    double smoothEnergySumMax;
    double energySumMax;
    
    if(grueneisen != 1)
    {
        for(int i=0; i<binDim[0]; i++)
        {
            for(int j=0; j<binDim[1]; j++)
            {
                for(int k=0; k<binDim[2]; k++)
                {
                    energySum[k+binDim[2]*(j+binDim[1]*i)] *= grueneisen; 
                }
            }
        }
    }
    
    // Obtain the max absolute value of 'energySum'
    if(restoreMagnitude)
    {
        double* eSumMax = std::max_element(energySum,energySum+(binDim[0]*binDim[1]*binDim[2]));
        double* eSumMin = std::min_element(energySum,energySum+(binDim[0]*binDim[1]*binDim[2]));
        if(std::abs(*eSumMin) > *eSumMax)
            energySumMax = std::abs(*eSumMin);
        else
            energySumMax = *eSumMax;
    }

    fftw_plan plan;
    fftw_complex in[binDim[0]*binDim[1]*binDim[2]];
    fftw_complex out[binDim[0]*binDim[1]*binDim[2]];

    // Fill 'in' with 'energySum' values, 
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            for(int k=0; k<binDim[2]; k++)
            {
                in[k+binDim[2]*(j+binDim[1]*i)][0] = energySum[k+binDim[2]*(j+binDim[1]*i)]; // Real
                in[k+binDim[2]*(j+binDim[1]*i)][1] = 0; // Imaginary
            }
        }
    }
        
    // Apply DFT from 'in' and write to 'out'
    plan = fftw_plan_dft_3d(binDim[0],binDim[1],binDim[2],in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
        
    // Which window filter type to multiply by ('smooth' in kWave uses rotated window by default)
    if( rotWindow )
        RotWindowFilter(binDim,out); // Rotation method to generate window filter
    else
        WindowFilter(binDim,out); // Multiply 'out' by window filter
    
        
    // Inverse DFT transformation and overwrite 'in'
    // There will be some small difference due to floating point precision when compared to matlab
    plan = fftw_plan_dft_3d(binDim[0],binDim[1],binDim[2],out,in,FFTW_BACKWARD,FFTW_ESTIMATE); 
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // Fill overwrite 'energySum' (a double array) with 'in' (an fftw_complex), also normalise
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            for(int k=0; k<binDim[2]; k++)
            {
                int index = k+binDim[2]*(j+binDim[1]*i);
                int indexTr = i+binDim[0]*(j+binDim[1]*k);             // We need to transpose data
                                                                      // into order K-Wave expects
                energySum[indexTr] = in[index][0];
                energySum[indexTr] /= (binDim[0]*binDim[1]*binDim[2]); // Normalise transformation
            }
        }
    }
    
    // Restore magnitude to before fftw transform
    if(restoreMagnitude)
    {
        // Get maximum absolute value of smoothed energySum
        double* smoothESumMax = std::max_element(energySum,energySum+(binDim[0]*binDim[1]*binDim[2]));
        double* smoothESumMin = std::min_element(energySum,energySum+(binDim[0]*binDim[1]*binDim[2])); 
        if(std::abs(*smoothESumMin) > *smoothESumMax)
            smoothEnergySumMax = std::abs(*smoothESumMin);
        else
            smoothEnergySumMax = *smoothESumMax;

        double restoreMag = (energySumMax)/(smoothEnergySumMax);

        // Multiply by each element to restore amplitude
        for(int i=0; i<binDim[0]; i++)
        {
            for(int j=0; j<binDim[1]; j++)
            {
                for(int k=0; k<binDim[2]; k++)
                {
                    int indexTr = i+binDim[0]*(j+binDim[1]*k);
                    energySum[indexTr] *= restoreMag; // Restore magnitude
                }
            }
        }
    }
    
    fftw_cleanup();

    // Convert data from double to float as per k-Wave input    
    float* fenergySum = new float[binDim[0]*binDim[1]*binDim[2]]{0}; // Float of data to be written
    std::copy(energySum, energySum+(binDim[0]*binDim[1]*binDim[2]), fenergySum);
    
    return fenergySum;
}


void KWaveInput::WriteHeader()
{
    time_t ttime = time(0);
    char* dt = ctime(&ttime);
    
    const char *createdBy[2] = { "created_by", "k-Wave 1.3" };
    //const char *creationDate[2] = { "creation_date", "31-Jan-2022-17-21-07" };
    const char *creationDate[2] = { "creation_date", dt };
    const char *fileDescription[2] = { "file_description", "Input data created for LhARA via KWaveInput." };
    const char *fileType[2] = { "file_type", "input" };
    const char *majorVersion[2] = { "major_version", "1" };
    const char *minorVersion[2] = { "minor_version", "2" };
    
    H5LTset_attribute_string(fileId, "/", createdBy[0], createdBy[1]);
    H5LTset_attribute_string(fileId, "/", creationDate[0], creationDate[1]);
    H5LTset_attribute_string(fileId, "/", fileDescription[0], fileDescription[1]);
    H5LTset_attribute_string(fileId, "/", fileType[0], fileType[1]);
    H5LTset_attribute_string(fileId, "/", majorVersion[0], majorVersion[1]);
    H5LTset_attribute_string(fileId, "/", minorVersion[0], minorVersion[1]);
}

// For calculating sensor indices (assuming matlab indexing for input which starts from 1)
// Implementation of sensor.mask() and find(sensor.mask) in matlab
std::vector<unsigned long> KWaveInput::SensorMaskIndices(unsigned long nx, int senMinX, int senMaxX, int senDx,
                                                         unsigned long ny, int senMinY, int senMaxY, int senDy,
                                                         unsigned long nz, int senMinZ, int senMaxZ, int senDz)
{
    std::vector<unsigned long> sensorIndices;
    unsigned long index;
    
    // Bad way to handle when senMin = senMax
    if(senDx == 0) 
        senDx = 1;
    if(senDy == 0)
        senDy = 1;
    if(senDz == 0)
        senDz = 1;
    
    // Calculate the index for sensors
    for(int senZ = senMinZ; senZ <= senMaxZ; senZ += senDz)
    {
        for(int senY = senMinY; senY <= senMaxY; senY += senDy)
        {
            for(int senX = senMinX; senX <= senMaxX; senX += senDx)
            {
                //std::cout << senZ << std::endl;
                index = ( nx * ny * ( senZ - 1 ) ) + ( nx * ( senY - 1 ) + senX );
                sensorIndices.push_back(index);
                //std::cout << index << std::endl;
            }
        }
    }
    
    return sensorIndices;
}

// Calculate dt
float KWaveInput::CalculateDT(float cfl, float* gridSpacing, float speed)
{
    int dim = sizeof(gridSpacing)/sizeof(gridSpacing[0]) + 1;   // calculate number of dimensions, 
                                                                // though for now code only works for 3d
    float minGrid = -1;
    
    for(int i=0; i<dim; i++)
    {
        if(minGrid == -1)
            minGrid = gridSpacing[i];
        else if(minGrid > gridSpacing[i])
            minGrid = gridSpacing[i];
    }
    float dt = (cfl * minGrid) / speed;
    return dt;
}

// Calculate Nt
unsigned long KWaveInput::CalculateNT(float tend, float dt)
{
    unsigned long Nt = floor(tend / dt) + 1;
    return Nt;
}


// For writing k-wave 1D float variable
void KWaveInput::Write1DFloat(std::string datasetName, float val)
{
    float value = val;
    
    // Attributes
    const char *dataType[2] = { "data_type", "float" };
    const char *domainType[2] = { "domain_type", "real" };

    hsize_t dimsf[3] = {1,1,1}; // 1D variable
    hid_t dataspace = H5Screate_simple(3, dimsf, nullptr);

    hid_t chunkProperty = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[3] = {1,1,1};
    herr_t status = H5Pset_chunk(chunkProperty, 3, chunk);

    const hid_t datasetType = H5T_NATIVE_FLOAT; // 32 bit floating point
    const hid_t pGroup = H5P_DEFAULT;

    hid_t dataset = H5Dcreate(fileId,datasetName.c_str(),datasetType,dataspace,H5P_DEFAULT,chunkProperty,H5P_DEFAULT);
    H5Dwrite(dataset, datasetType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    
    H5LTset_attribute_string(fileId,datasetName.c_str(),dataType[0],dataType[1]);
    H5LTset_attribute_string(fileId,datasetName.c_str(),domainType[0],domainType[1]);
}

// For writing k-wave 1D long variable
void KWaveInput::Write1DLong(std::string datasetName, unsigned long val)
{
    unsigned long value = val;
    const char *dataType[2] = { "data_type", "long" };
    const char *domainType[2] = { "domain_type", "real" };

    hsize_t dimsf[3] = {1,1,1}; // 1D variable
    hid_t dataspace = H5Screate_simple(3, dimsf, nullptr);

    hid_t chunkProperty = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[3] = {1,1,1};
    herr_t status = H5Pset_chunk(chunkProperty, 3, chunk);

    const hid_t datasetType = H5T_STD_U64LE; // 64 bit unsigned integer
    const hid_t pGroup = H5P_DEFAULT;

    hid_t dataset = H5Dcreate(fileId,datasetName.c_str(),datasetType,dataspace,H5P_DEFAULT,chunkProperty,H5P_DEFAULT);
    H5Dwrite(dataset, datasetType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &value);
    
    H5LTset_attribute_string(fileId,datasetName.c_str(),dataType[0],dataType[1]);
    H5LTset_attribute_string(fileId,datasetName.c_str(),domainType[0],domainType[1]);
}

// For writing k-wave mult dim long variable (i.e. sensor_mask_index)
void KWaveInput::WriteMultDLong(std::string datasetName, std::vector<unsigned long>& indicesVec)
{
    const char *dataType[2] = { "data_type", "long" };
    const char *domainType[2] = { "domain_type", "real" };

    hsize_t rank = indicesVec.size(); // Number of indices
    
    // H5DWrite doesn't accept vectors so convert to array
    unsigned long indicesArr[rank];
    for( int i=0; i<rank; i++)
        indicesArr[i] = indicesVec[i];
    
    hsize_t dimsf[3] = {1,1,rank}; // 1D variable
    hid_t dataspace = H5Screate_simple(3, dimsf, nullptr);

    hid_t chunkProperty = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[3] = {1,1,rank};
    herr_t status = H5Pset_chunk(chunkProperty, 3, chunk);

    const hid_t datasetType = H5T_STD_U64LE; // 64 bit unsigned integer
    const hid_t pGroup = H5P_DEFAULT;

    hid_t dataset = H5Dcreate(fileId,datasetName.c_str(),datasetType,dataspace,H5P_DEFAULT,chunkProperty,H5P_DEFAULT);
    
    H5Dwrite(dataset, datasetType, H5S_ALL, H5S_ALL, H5P_DEFAULT, &indicesArr);
    
    H5LTset_attribute_string(fileId,datasetName.c_str(),dataType[0],dataType[1]);
    H5LTset_attribute_string(fileId,datasetName.c_str(),domainType[0],domainType[1]);
}

void KWaveInput::WritePData(std::string datasetName, int nx, int ny, int nz, float* data)
{
    const char *dataType[2] = { "data_type", "float" };
    const char *domainType[2]= { "domain_type", "real" };

    hsize_t dimsf[3];
    dimsf[0] = nx;
    dimsf[1] = ny;
    dimsf[2] = nz;
    hid_t dataspace = H5Screate_simple(3, dimsf, nullptr);
    
    hid_t chunkProperty = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[3];
    chunk[0] = 1;
    chunk[1] = ny;
    chunk[2] = nz;
    herr_t status = H5Pset_chunk(chunkProperty, 3, chunk);

    const hid_t datasetType = H5T_NATIVE_FLOAT; // 32 bit floating point
    const hid_t pGroup = H5P_DEFAULT;

    hid_t dataset = H5Dcreate(fileId,datasetName.c_str(),datasetType,dataspace,H5P_DEFAULT,chunkProperty,H5P_DEFAULT);
    H5Dwrite(dataset, datasetType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
    H5LTset_attribute_string(fileId,datasetName.c_str(),dataType[0],dataType[1]);
    H5LTset_attribute_string(fileId,datasetName.c_str(),domainType[0],domainType[1]);
}

void KWaveInput::WriteFile()
{
    std::cout << "Beginning to write K-Wave Input File." << std::endl;
    
    if( enableSmooth )
        fenergyData = Smooth(energyData, binDim, grueneisen, restoreMagnitude, rotWindow);
    else
    {
        fenergyData = new float[binDim[0]*binDim[1]*binDim[2]]{0};
        std::copy(energyData, energyData+(binDim[0]*binDim[1]*binDim[2]), fenergyData);
    }
    
    float gridSpacing[3];
    gridSpacing[0] = dx;
    gridSpacing[1] = dy;
    gridSpacing[2] = dz;
    
    float dt = CalculateDT(cfl, gridSpacing, mediumSoundSpeed);
    unsigned long Nt = CalculateNT(tEnd, dt);

    // Writing Variables to Input File
    Write1DLong("Nt",Nt);
    Write1DLong("Nx",Nx);
    Write1DLong("Ny",Ny);
    Write1DLong("Nz",Nz);

    Write1DLong("absorbing_flag",0);
    Write1DLong("axisymmetric_flag",0);
    
    Write1DFloat("c0",mediumSoundSpeed);
    Write1DFloat("c_ref",cRef);

    Write1DFloat("dt",dt);
    Write1DFloat("dx",dx);
    Write1DFloat("dy",dy);
    Write1DFloat("dz",dz);

    Write1DLong("elastic_flag",0);
    Write1DLong("nonlinear_flag",0);
    Write1DLong("nonuniform_grid_flag",0);

    Write1DLong("p0_source_flag",1);
    WritePData("p0_source_input", Nz, Ny, Nx, fenergyData);
    Write1DLong("p_source_flag",0);

    Write1DFloat("pml_x_alpha",pmlAlphaX);
    Write1DLong("pml_x_size",pmlSizeX);
    Write1DFloat("pml_y_alpha",pmlAlphaY);
    Write1DLong("pml_y_size",pmlSizeY);
    Write1DFloat("pml_z_alpha",pmlAlphaZ);
    Write1DLong("pml_z_size",pmlSizeZ);
    
    Write1DFloat("rho0",1);
    Write1DFloat("rho0_sgx",1);
    Write1DFloat("rho0_sgy",1);
    Write1DFloat("rho0_sgz",1);

    WriteMultDLong("sensor_mask_index",sensorMask);
    
    Write1DLong("sensor_mask_type",0);
    Write1DLong("sxx_source_flag",0);
    Write1DLong("sxy_source_flag",0);
    Write1DLong("sxz_source_flag",0);
    Write1DLong("syy_source_flag",0);
    Write1DLong("syz_source_flag",0);
    Write1DLong("szz_source_flag",0);

    Write1DLong("transducer_source_flag",0);

    Write1DLong("ux_source_flag",0);
    Write1DLong("uy_source_flag",0);
    Write1DLong("uz_source_flag",0);
}
