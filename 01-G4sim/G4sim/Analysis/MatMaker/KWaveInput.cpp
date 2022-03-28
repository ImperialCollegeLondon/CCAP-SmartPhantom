#include "KWaveInput.h"

KWaveInput::KWaveInput()
{}

KWaveInput::KWaveInput(const std::string& fn)
{
    grueneisen = 1.;    // Set default value 
    fileId = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    WriteHeader();
}

KWaveInput::~KWaveInput()
{
    std::cout << "Beginning destructor for KWaveInput." << std::endl;
    delete fenergyDensityData;
    H5Fclose(fileId);
    std::cout << "End of destructor for KWaveInput." << std::endl;    
}


void KWaveInput::circshift(double *out, const double *in, int xdim, int ydim, int xshift, int yshift)
{
    /*
        Shift array circularly, taken from
        https://stackoverflow.com/questions/5915125/fftshift-ifftshift-c-c-source-code
        
        Input:      out    -- output array
                    in     -- input array
                    xdim   -- dimension size of x
                    ydim   -- dimension size of y
                    xshift -- shift in x 
                    yshift -- shift in y
    */

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

std::vector<double> KWaveInput::Linspace(double first, double last, int len)
{
    /*
        Return an evenly spaced 1-d grid of doubles taken from 
        https://rncarpio.github.io/linterp/
        
        Input:      first -- Initial value
                    last  -- Final value
                    len   -- Number of values
                    
        Return:     Evenly spaced 1-d grid
    */

    std::vector<double> result(len);
    double step = (last - first) / (len - 1);
    for(int i=0; i<len; i++)
    {
        result[i] = first + i*step;
    }
    return result;
}

double* KWaveInput::Blackman(int wSize)
{
    /*
        Create Blackman window. Meant to produce the same result as 
        blackman from scipy.signal.
        
        Input:      wSize -- Number of points in output window
                    
        Return:     Blackman window
    */

    int windowSize = wSize;
    window = new double[windowSize];
    
    if( windowSize%2 == 0) // If even, return periodic window
        windowSize += 1;
    
    double stepsize = (double)windowSize/(windowSize - 1);
    const double pi = std::acos(-1);
    
    // Default parameters
    double a0 = 0.42;
    double a1 = 0.5;
    double a2 = 0.08;
    
    double n = 0;
    double maxVal = 0;
    
    for(int i=0; i<wSize; i++) // counter is wrt wSize to handle even and odd as specified in smooth.m from k-wave
    {
        double blackman = a0 - a1*std::cos(2*pi*n/windowSize) + a2*std::cos(4*pi*n/windowSize);
        window[i] = blackman;
        n += stepsize;
        if(maxVal < blackman)
            maxVal = blackman;
    }

    return window;
}

void KWaveInput::CreateFilter(int* adjustData, double* filterCombined)
{
    /*
        Create Blackman window filter by combining the dimensions
        
        Input:      adjustData      -- Adjusted voxel dimensions 
                    filterCombined  -- array of filter                    
    */

    if(binDimSize == 2)
    {
        double filterX[adjustData[0]] = {0};
        double filterY[adjustData[1]] = {0};

        // Getting filter values for each dimension separately
        // **************************************************************
        for(int axis=0; axis<binDimSize; axis++)
        {
            double* blackman = Blackman(adjustData[axis]);
            for(int i=0; i<adjustData[axis]; i++)
            {
                if(axis==0)
                    filterX[i] = blackman[i];
                else if(axis==1)
                    filterY[i] = blackman[i];
            }
            delete window;
        }

        // Multiply filters together
        // **************************************************************
        for(int i=0; i<adjustData[0]; i++)
        {
            for(int j=0; j<adjustData[1]; j++)
            {
                int index = j+adjustData[1]*i;
                filterCombined[index] = std::abs(filterX[i]*filterY[j]);
            }
        }
    }
    else if(binDimSize == 3)
    {
        double filterX[adjustData[0]] = {0};
        double filterY[adjustData[1]] = {0};
        double filterZ[adjustData[2]] = {0};
        
        // Getting filter values for each dimension separately
        // **************************************************************
        for(int axis=0; axis<binDimSize; axis++)
        {
            double* blackman = Blackman(adjustData[axis]);
            for(int i=0; i<adjustData[axis]; i++)
            {
                if(axis==0)
                    filterX[i] = blackman[i];
                else if(axis==1)
                    filterY[i] = blackman[i];
                else if(axis==2)
                    filterZ[i] = blackman[i];
            }
            delete window;
        }
        
        // Multiply filters together
        // **************************************************************
        for(int i=0; i<adjustData[0]; i++)
        {
            for(int j=0; j<adjustData[1]; j++)
            {
                for(int k=0; k<adjustData[2]; k++)
                {
                    int index = k + adjustData[2]*(j+adjustData[1]*i);
                    filterCombined[index] = std::abs(filterX[i]*filterY[j]*filterZ[k]);
                }
            }
        }
    }
    else
    {
        std::cout << "Unexpected dimension size of " << binDimSize << ", provided in KWaveInput::CreateFilter()." << std::endl;
    }
}

void KWaveInput::SliceFilter(int* nArr, int* adjustData, double* filterCombined, double* filterSliced)
{
    /*
        Slice the filter size to match voxel size (nArr), quite messy, needs a better implementation.
        
        Input:      nArr            -- Voxel Dimensions
                    adjustData      -- Adjusted voxel dimensions 
                    filterCombined  -- array of filter
                    filterSliced    -- array to hold sliced filter using nArr dimensions                    
    */
    
    int sliceIndex = 0;

    if(binDimSize == 2)
    {
        for(int i=0; i<nArr[0]; i++)
        {
            for(int j=0; j<nArr[1]; j++)
            {
                int index = j+nArr[1]*i;
                filterSliced[index] = filterCombined[sliceIndex];
                sliceIndex++;
            }
            if(nArr[1]!=adjustData[1])
                sliceIndex++;
        }
    }
    else if(binDimSize == 3)
    {
        for(int i=0; i<nArr[0]; i++)
        {
            for(int j=0; j<nArr[1]; j++)
            {
                for(int k=0; k<nArr[2]; k++)
                {
                    int index = k + nArr[2]*(j+nArr[1]*i);
                    filterSliced[index] = filterCombined[sliceIndex];
                    sliceIndex++;
                }
                if(nArr[2]!=adjustData[2])
                    sliceIndex++;
            }
            if(nArr[1]!=adjustData[1])
                sliceIndex += adjustData[2];
        }
    }
    else
    {
        std::cout << "Unexpected dimension size of " << binDimSize << ", provided in KWaveInput::SliceFilter()." << std::endl;
    }
}

std::vector<double> KWaveInput::CreateRadiusGrid(int* adjustData, double& radius)
{
    /*
        Create radius grid Blackman window filter and applied to 3D energy data
        
        Input:      adjustData -- Adjusted voxel dimensions
                    radius     -- Radius
                    
        Return:     Filled radius grid
    */
    
    std::vector<double> radiusGrid;
    
    if(binDimSize == 2)
    {
        radiusGrid.reserve(adjustData[0]*adjustData[1]);
        
        std::vector<double> xx = Linspace(-radius,radius,adjustData[0]);
        std::vector<double> yy = Linspace(-radius,radius,adjustData[1]);

        // Set up grid based on radius
        for(int i=0; i<adjustData[0]; i++)
        {
            for(int j=0; j<adjustData[1]; j++)
            {
                int index = j+adjustData[1]*i;
                double rr = std::sqrt( xx[i]*xx[i] + yy[j]*yy[j] );
                if( rr > radius)
                    rr = radius;
                radiusGrid[index] = rr;
            }
        }
    }
    else if(binDimSize == 3)
    {
        radiusGrid.reserve(adjustData[0]*adjustData[1]*adjustData[2]);
        
        std::vector<double> xx = Linspace(-radius,radius,adjustData[0]);
        std::vector<double> yy = Linspace(-radius,radius,adjustData[1]);
        std::vector<double> zz = Linspace(-radius,radius,adjustData[2]);

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
                    radiusGrid[index] = rr;
                }
            }
        }
    }
    else
    {
        std::cout << "Unexpected dimension size of " << binDimSize << ", provided in KWaveInput::CreateRadiusGrid()." << std::endl;
    }
    
    return radiusGrid;
}

std::vector<double> KWaveInput::CreateWindow(int& L)
{
    /*
        Create linear Blackman window
        
        Input:      L -- Number of elements in window
                    
        Return:     Blackman window as a vector
    */
    
    std::vector<double> winLin;
    winLin.reserve(L);
    double* bWin = Blackman(L);

    for(int i=0; i<L; i++)
        winLin[i] = bWin[i];

    return winLin;
}

void KWaveInput::CreateRotFilter(int* nArr, int* adjustData, int& L, double& radius, 
                                    std::vector<double>& radiusGrid, std::vector<double>& winLin,
                                    double* rotWin)
{
    /*
        Create rotation window filter from interpolated values
        
        Input:      nArr       -- Voxel dimensions
                    adjustData -- Adjusted voxel dimensions
                    L          -- Number of elements in window
                    radius     -- Radius
                    radiusGrid -- Radius Grid
                    winLin     -- Blackman window of size L
                    rotWin     -- Container for rotation window filter size of nArr                    
    */

    std::vector<double> ll = Linspace(-radius,radius,L);
    
    // Set up interpolator parameters
    // **************************************************************
    std::vector<std::vector<double>::iterator> gridIterList;
    gridIterList.push_back(ll.begin());
    array<int,1> gridSizes;
    gridSizes[0] = L;
    InterpMultilinear<1, double> 
    interp_ML(gridIterList.begin(),gridSizes.begin(),winLin.data(),winLin.data()+L);

    // Write interpolated values to rotation window
    // **************************************************************
    if( binDimSize == 2)
    {
        int rIndex = 0;

        for(int i=0; i<nArr[0]; i++)
        {
            for(int j=0; j<nArr[1]; j++)
            {
                int rotIndex = j+nArr[1]*i;
                array<double,1> args = {radiusGrid[rIndex]};
                rotWin[rotIndex] = std::abs(interp_ML.interp(args.begin()));
                rIndex++;
            }
            if(nArr[1]!=adjustData[1])
                rIndex++;
        }
    }
    else if( binDimSize == 3)
    {
        int rIndex = 0;

        for(int i=0; i<nArr[0]; i++)
        {
            for(int j=0; j<nArr[1]; j++)
            {
                for(int k=0; k<nArr[2]; k++)
                {
                    int rotIndex = k + nArr[2]*(j+nArr[1]*i);
                    array<double,1> args = {radiusGrid[rIndex]};
                    rotWin[rotIndex] = std::abs(interp_ML.interp(args.begin()));
                    rIndex++;
                }
                if(nArr[2]!=adjustData[2])
                    rIndex++;
            }
            if(nArr[1]!=adjustData[1])
                rIndex += adjustData[2];
        }   
    }
    else
    {
        std::cout << "Unexpected dimension size of " << binDimSize << ", provided in KWaveInput::CreateRotFilter()." << std::endl;
    }
}

void KWaveInput::WindowFilter(int* nArr, fftw_complex* data)
{
    /*
        Create a Blackman window filter and applied to 3D energy data
        Can give different results compared to Matlab
        
        Input:      nArr -- Voxel Dimensions
                    data -- FFT energy data
    */
    
    // Check dimensions are correct
    // **************************************************************
    if(binDimSize != 3)
    {
        std::cout << "KWaveInput::WindowFilter() only accepts 3D, instead received a size = " << binDimSize << std::endl;
        return;
    }

    // To create Non-symmetrical window for even dimensions
    // **************************************************************
    int adjustData[binDimSize];
    
    for(int i=0; i<binDimSize; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }

    // Functions to create filter, slice to appropriate dimensions, and apply ifftshift
    // **************************************************************
    double filterCombined[adjustData[0]*adjustData[1]*adjustData[2]];
    double filterSliced[nArr[0]*nArr[1]*nArr[2]];
    double ifftFilterCombined[nArr[0]*nArr[1]*nArr[2]];

    CreateFilter(adjustData, filterCombined);
    SliceFilter(nArr, adjustData, filterCombined, filterSliced);
    ifftshift(ifftFilterCombined, filterSliced, nArr[0], nArr[1], nArr[2]);    
    
    // Apply the filter to the data 
    // **************************************************************
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int index = k+nArr[2]*(j+nArr[1]*i);
                data[index][0] *= ifftFilterCombined[index]; // Real
                data[index][1] *= ifftFilterCombined[index]; // Imaginary
            }
        }
    }
}

void KWaveInput::WindowFilterAS(int* nArr, fftw_complex* data)
{
    /*
        Create a Blackman window filter and applied to 2D energy data (axisymmetric)
        Can give different results compared to Matlab
        
        Input:      nArr -- Voxel Dimensions
                    data -- FFT energy data
    */

    // To create Non-symmetrical window for even dimensions
    // **************************************************************
    if(binDimSize != 2)
    {
        std::cout << "KWaveInput::WindowFilterAS() only accepts 2D, instead received a size = " << binDimSize << std::endl;
        return;
    }

    // To create Non-symmetrical window for even dimensions
    // **************************************************************
    int adjustData[binDimSize];
    
    for(int i=0; i<binDimSize; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }

    // Functions to create filter, slice to appropriate dimensions, and apply ifftshift
    // **************************************************************
    double filterCombined[adjustData[0]*adjustData[1]];
    double filterSliced[nArr[0]*nArr[1]];
    double ifftFilterCombined[nArr[0]*nArr[1]];

    CreateFilter(adjustData, filterCombined);
    SliceFilter(nArr, adjustData, filterCombined, filterSliced);
    ifftshift(ifftFilterCombined, filterSliced, nArr[0], nArr[1]);    
    
    // Apply the filter to the data 
    // **************************************************************
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            int index = j+nArr[1]*i;
            data[index][0] *= ifftFilterCombined[index]; // Real
            data[index][1] *= ifftFilterCombined[index]; // Imaginary
        }
    }
}

void KWaveInput::RotWindowFilter(int* nArr, fftw_complex* data)
{
    /*
        Create a Blackman window filter via rotation method and applied to 3D energy data
        
        Input:      nArr -- Voxel Dimensions
                    data -- FFT energy data
    */

    // Check dimensions are correct
    // **************************************************************
    if(binDimSize != 3)
    {
        std::cout << "KWaveInput::RotWindowFilter() only accepts 3D, instead received a size = " << binDimSize << std::endl;
        return;
    }
    
    // Adjust dimensions to account for nonsymmetric filter if dimension has even number of terms
    // **************************************************************
    int adjustData[binDimSize];
    for(int i=0; i<binDimSize; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }
    
    int L = *std::max_element(adjustData,adjustData+binDimSize);
    double radius = (L-1)/2.;
    
    std::vector<double> radiusGrid = CreateRadiusGrid(adjustData, radius);
    std::vector<double>     winLin = CreateWindow(L);
        
    double rotWin[nArr[0]*nArr[1]*nArr[2]] = {0};
    double ifftRotWin[nArr[0]*nArr[1]*nArr[2]] = {0};

    // Functions to create rotation filter, slice to appropriate dimensions, and apply ifftshift
    // **************************************************************
    CreateRotFilter(nArr, adjustData, L, radius, radiusGrid, winLin, rotWin);
    ifftshift(ifftRotWin, rotWin, nArr[0], nArr[1], nArr[2]);    
    
    // Apply the filter to the data 
    // **************************************************************
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            for(int k=0; k<nArr[2]; k++)
            {
                int index = k+nArr[2]*(j+nArr[1]*i);
                data[index][0] *= ifftRotWin[index]; // Real
                data[index][1] *= ifftRotWin[index]; // Imaginary
            }
        }
    }
}

void KWaveInput::RotWindowFilterAS(int* nArr, fftw_complex* data)
{    
    /*
        Create a Blackman window filter via rotation method and applied to 2D energy data (axisymmetric)
        
        Input:      nArr -- Voxel Dimensions
                    data -- FFT energy data
    */
    
    // Check dimensions are correct
    // **************************************************************
    if(binDimSize != 2)
    {
        std::cout << "KWaveInput::RotWindowFilterAS() only accepts 2D, instead received a size = " << binDimSize << std::endl;
        return;
    }

    // Adjust dimensions to account for nonsymmetric filter if dimension has even number of terms
    // **************************************************************
    int adjustData[binDimSize];
    for(int i=0; i<binDimSize; i++)
    {
        if(nArr[i]%2 == 0)
            adjustData[i] = nArr[i] + 1; // Add one to dimension if even to make nonsymmetric filter
        else
            adjustData[i] = nArr[i];
    }

    int L = *std::max_element(adjustData,adjustData+binDimSize);
    double radius = (L-1)/2.;

    std::vector<double> radiusGrid = CreateRadiusGrid(adjustData, radius);
    std::vector<double>     winLin = CreateWindow(L);
        
    double rotWin[nArr[0]*nArr[1]] = {0};
    double ifftRotWin[nArr[0]*nArr[1]] = {0};

    // Functions to create rotation filter, slice to appropriate dimensions, and apply ifftshift
    // **************************************************************
    CreateRotFilter(nArr, adjustData, L, radius, radiusGrid, winLin, rotWin);
    ifftshift(ifftRotWin, rotWin, nArr[0], nArr[1]);    
    
    // Apply the filter to the data 
    // **************************************************************
    for(int i=0; i<nArr[0]; i++)
    {
        for(int j=0; j<nArr[1]; j++)
        {
            int index = j+nArr[1]*i;
            data[index][0] *= ifftRotWin[index]; // Real
            data[index][1] *= ifftRotWin[index]; // Imaginary
        }
    }
}

float* KWaveInput::Smooth(double* energyDensitySum, int* binDim, bool restoreMagnitude, bool rotWindow)
{
    /*
        Apply smooth procedure to 3D energy data
        
        Input:      energyDensitySum -- Binned energy data
                    binDim           -- Voxel dimensions
                    restoreMagnitude -- Boolean on whether to restore magnitude to transformed data
                    rotWindow        -- Whether to use the rotational method to generate window filter
                    
        Return:     Smoothed energy data
    */

    double smoothEnergyDensitySumMax;
    double energyDensitySumMax;

    // Obtain the max absolute value of 'energyDensitySum'
    // **************************************************************
    if(restoreMagnitude)
    {
        double* eSumMax = std::max_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]*binDim[2]));
        double* eSumMin = std::min_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]*binDim[2]));
        if(std::abs(*eSumMin) > *eSumMax)
            energyDensitySumMax = std::abs(*eSumMin);
        else
            energyDensitySumMax = *eSumMax;
    }

    fftw_plan plan;
    fftw_complex in[binDim[0]*binDim[1]*binDim[2]];
    fftw_complex out[binDim[0]*binDim[1]*binDim[2]];

    // Fill 'in' with 'energyDensitySum' values, 
    // **************************************************************
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            for(int k=0; k<binDim[2]; k++)
            {
                int index = k+binDim[2]*(j+binDim[1]*i);
                in[index][0] = energyDensitySum[index];            // Real
                in[index][1] = 0;                           // Imaginary
            }
        }
    }
        
    // Apply DFT from 'in' and write to 'out'
    // **************************************************************
    plan = fftw_plan_dft_3d(binDim[0],binDim[1],binDim[2],in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
        
    // Which window filter type to multiply by ('smooth' in kWave uses rotated window by default)
    // **************************************************************
    if( rotWindow )
        RotWindowFilter(binDim,out);    // Rotation method to generate window filter
    else
        WindowFilter(binDim,out);       // Multiply 'out' by window filter
            
    // Inverse DFT transformation and overwrite 'in'
    // **************************************************************
    plan = fftw_plan_dft_3d(binDim[0],binDim[1],binDim[2],out,in,FFTW_BACKWARD,FFTW_ESTIMATE); 
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    
    // Fill overwrite 'energyDensitySum' (a double array) with 'in' (an fftw_complex), also normalise
    // **************************************************************
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            for(int k=0; k<binDim[2]; k++)
            {
                int index = k+binDim[2]*(j+binDim[1]*i);
                int indexTr = i+binDim[0]*(j+binDim[1]*k);            // We need to transpose data
                                                                      // into order K-Wave expects
                energyDensitySum[indexTr] = in[index][0];
                energyDensitySum[indexTr] /= (binDim[0]*binDim[1]*binDim[2]); // Normalise transformation
            }
        }
    }
    
    // Restore magnitude to before fourier transform
    // **************************************************************
    if(restoreMagnitude)
    {
        // Get maximum absolute value of smoothed energyDensitySum
        // **************************************************************
        double* smoothESumMax = std::max_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]*binDim[2]));
        double* smoothESumMin = std::min_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]*binDim[2])); 
        if(std::abs(*smoothESumMin) > *smoothESumMax)
            smoothEnergyDensitySumMax = std::abs(*smoothESumMin);
        else
            smoothEnergyDensitySumMax = *smoothESumMax;

        double restoreMag = (energyDensitySumMax)/(smoothEnergyDensitySumMax);

        // Multiply by each element to restore amplitude
        // **************************************************************
        for(int i=0; i<binDim[0]; i++)
        {
            for(int j=0; j<binDim[1]; j++)
            {
                for(int k=0; k<binDim[2]; k++)
                {
                    int indexTr = i+binDim[0]*(j+binDim[1]*k);
                    energyDensitySum[indexTr] *= restoreMag; // Restore magnitude
                }
            }
        }
    }
    
    fftw_cleanup();

    // Convert data from double to float as per k-Wave input    
    // **************************************************************
    float* fenergyDensitySum = new float[binDim[0]*binDim[1]*binDim[2]]{0}; // Float of data to be written
    std::copy(energyDensitySum, energyDensitySum+(binDim[0]*binDim[1]*binDim[2]), fenergyDensitySum);
    
    return fenergyDensitySum;
}

float* KWaveInput::SmoothAS(double* energyDensitySum, int* binDim, bool restoreMagnitude, bool rotWindow)
{
    /*
        Apply smooth procedure to 2D energy density data (axisymmetric)
        
        Input:      energyDensitySum -- Binned energy density data
                    binDim           -- Voxel dimensions
                    restoreMagnitude -- Boolean on whether to restore magnitude to transformed data
                    rotWindow        -- Whether to use the rotational method to generate window filter
                    
        Return:     Smoothed energy density data
    */

    double smoothEnergyDensitySumMax;
    double energyDensitySumMax;

    // Obtain the max absolute value of 'energyDensitySum'
    // **************************************************************
    if(restoreMagnitude)
    {
        double* eSumMax = std::max_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]));
        double* eSumMin = std::min_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]));
        if(std::abs(*eSumMin) > *eSumMax)
            energyDensitySumMax = std::abs(*eSumMin);
        else
            energyDensitySumMax = *eSumMax;
    }
    
    fftw_plan plan;
    fftw_complex in[binDim[0]*binDim[1]];
    fftw_complex out[binDim[0]*binDim[1]];

    // Fill 'in' with 'energyDensitySum' values, 
    // **************************************************************
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            int index = j + binDim[1] * i;
            in[index][0] = energyDensitySum[index];    // Real
            in[index][1] = 0;                          // Imaginary
        }
    }
            
    // Apply DFT from 'in' and write to 'out'
    // **************************************************************
    plan = fftw_plan_dft_2d(binDim[0],binDim[1],in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
        
    // Which window filter type to multiply by ('smooth' in kWave uses rotated window by default)
    // **************************************************************
    if( rotWindow )
        RotWindowFilterAS(binDim,out);  // Rotation method to generate window filter
    else
        WindowFilterAS(binDim,out);     // Multiply 'out' by window filter
            
    // Inverse DFT transformation and overwrite 'in'
    // **************************************************************
    plan = fftw_plan_dft_2d(binDim[0],binDim[1],out,in,FFTW_BACKWARD,FFTW_ESTIMATE); 
    fftw_execute(plan);
    fftw_destroy_plan(plan);
        
    // Fill overwrite 'energyDensitySum' (a double array) with 'in' (an fftw_complex), also normalise
    // **************************************************************
    for(int i=0; i<binDim[0]; i++)
    {
        for(int j=0; j<binDim[1]; j++)
        {
            int index = j+binDim[1]*i;
            energyDensitySum[index] = in[index][0];
            energyDensitySum[index] /= (binDim[0]*binDim[1]); // Normalise transformation
        }
    }
        
    // Restore magnitude to before fftw transform
    // **************************************************************
    if(restoreMagnitude)
    {
        // Get maximum absolute value of smoothed energyDensitySum
        // **************************************************************
        double* smoothESumMax = std::max_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1]));
        double* smoothESumMin = std::min_element(energyDensitySum,energyDensitySum+(binDim[0]*binDim[1])); 
    
        if(std::abs(*smoothESumMin) > *smoothESumMax)
            smoothEnergyDensitySumMax = std::abs(*smoothESumMin);
        else
            smoothEnergyDensitySumMax = *smoothESumMax;

        double restoreMag = (energyDensitySumMax)/(smoothEnergyDensitySumMax);

        // Multiply by each element to restore amplitude
        // **************************************************************
        for(int i=0; i<binDim[0]; i++)
        {
            for(int j=0; j<binDim[1]; j++)
            {
                int index = j+binDim[1]*i;
                energyDensitySum[index] *= restoreMag;
            }
        }
    }
    
    fftw_cleanup();

    // Convert data from double to float as per k-Wave input   
    // **************************************************************
    float* fenergyDensitySum = new float[binDim[0]*binDim[1]]{0}; // Float of data to be written
    std::copy(energyDensitySum, energyDensitySum+(binDim[0]*binDim[1]), fenergyDensitySum);
    
    return fenergyDensitySum;
}

void KWaveInput::WriteHeader()
{
    /*
        Writes header part for k-Wave input file
    */

    time_t ttime = time(0);
    char* dt = ctime(&ttime);
    
    const char *createdBy[2] = { "created_by", "k-Wave 1.3" };
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

std::vector<unsigned long> KWaveInput::SensorMaskIndices(unsigned long nx, float senMinX, float senMaxX, int senDx,
                                                         unsigned long ny, float senMinY, float senMaxY, int senDy,
                                                         unsigned long nz, float senMinZ, float senMaxZ, int senDz)
{
    /*
        Indices of point sensors (assuming matlab indexing for input which starts from 1)
        Implementation of sensor.mask() and find(sensor.mask)
        
        Input:      nx -- Voxel-X, senMinX -- Minimum-X, senMaxX -- Maximum-X, senDx -- Step-X
                    ny -- Voxel-Y, senMinY -- Minimum-Y, senMaxY -- Maximum-Y, senDy -- Step-Y
                    nz -- Voxel-Z, senMinZ -- Minimum-Z, senMaxZ -- Maximum-Z, senDz -- Step-Z
                    
        Return:     Point sensor indices
    */

    std::vector<unsigned long> sensorIndices;
    unsigned long index;
    
    // Bad way to handle when senMin = senMax
    // **************************************************************
    if(senDx == 0) 
        senDx = 1;
    if(senDy == 0)
        senDy = 1;
    if(senDz == 0)
        senDz = 1;
    
    // Round min/max
    // **************************************************************
    float minX = std::round(senMinX);
    float maxX = std::round(senMaxX);
    float minY = std::round(senMinY);
    float maxY = std::round(senMaxY);
    float minZ = std::round(senMinZ);
    float maxZ = std::round(senMaxZ);
        
    // Calculate the index for sensors
    // **************************************************************
    for(int senZ = minZ; senZ <= maxZ; senZ += senDz)
    {
        for(int senY = minY; senY <= maxY; senY += senDy)
        {
            for(int senX = minX; senX <= maxX; senX += senDx)
            {
                index = ( nx * ny * ( senZ - 1 ) ) + ( nx * ( senY - 1 ) + senX );
                sensorIndices.push_back(index);
            }
        }
    }

    return sensorIndices;
}

std::vector<unsigned long> KWaveInput::SensorMaskIndices(unsigned long nx, float senMinX, float senMaxX, int senDx,
                                                         unsigned long ny, float senMinY, float senMaxY, int senDy)
{
    /*
        Indices of point sensors
        
        Input:      nx -- Voxel-X, senMinX -- Minimum-X, senMaxX -- Maximum-X, senDx -- Step-X
                    ny -- Voxel-Y, senMinY -- Minimum-Y, senMaxY -- Maximum-Y, senDy -- Step-Y
                    
        Return:     Point sensor indices
    */

    std::vector<unsigned long> sensorIndices;
    unsigned long index;
    
    // Bad way to handle when senMin = senMax
    // **************************************************************
    if(senDx == 0) 
        senDx = 1;
    if(senDy == 0)
        senDy = 1;

    // Round min/max to int values
    // **************************************************************
    float minX = std::round(senMinX);
    float maxX = std::round(senMaxX);
    float minY = std::round(senMinY);
    float maxY = std::round(senMaxY);

    // Calculate the index for sensors
    // **************************************************************
    for(int senY = minY; senY <= maxY; senY += senDy)
    {
        for(int senX = minX; senX <= maxX; senX += senDx)
        {
            index = ( nx * ( senY - 1 ) + senX );
            sensorIndices.push_back(index);
        }
    }
    
    return sensorIndices;
}

float KWaveInput::CalculateDT(float cfl, float* gridSpacing, float speed)
{
    /*
        Calculate dt
        
        Input:      cfl         -- Courant Number
                    gridSpacing -- Grid spacing [m]
                    speed       -- Speed of sound
                    
        Return:     Value of dt
    */

    float minGrid = -1;
    
    for(int i=0; i<binDimSize; i++)
    {
        if(minGrid == -1)
            minGrid = gridSpacing[i];
        else if(minGrid > gridSpacing[i])
            minGrid = gridSpacing[i];
    }
    
    float dt = (cfl * minGrid) / speed;
    return dt;
}

unsigned long KWaveInput::CalculateNT(float tend, float dt)
{
    /*
        Calculates Nt
        
        Input:      tend -- End time
                    dt   -- time step
                    
        Return:     Value of Nt
    */

    unsigned long Nt = floor(tend / dt) + 1;
    return Nt;
}

void KWaveInput::Write1DFloat(std::string datasetName, float val)
{
    /*
        For writing k-wave 1D float variable to HDF5 file
        
        Input:      datasetName -- Name of variable
                    value       -- Value of variable
    */

    float value = val;
    
    // Attributes
    // **************************************************************
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

void KWaveInput::Write1DLong(std::string datasetName, unsigned long val)
{
    /*
        For writing k-wave 1D long variable to HDF5 file
        
        Input:      datasetName -- Name of variable
                    value       -- Value of variable
    */

    unsigned long value = val;

    // Attributes
    // **************************************************************    
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

void KWaveInput::WriteMultDLong(std::string datasetName, std::vector<unsigned long>& indicesVec)
{
    /*
        For writing k-wave vector of long variable (i.e. sensor_mask_index) to HDF5 file      

        Input:      datasetName -- Name of variable
                    value       -- Vector of values
    */

    const char *dataType[2] = { "data_type", "long" };
    const char *domainType[2] = { "domain_type", "real" };

    hsize_t rank = indicesVec.size(); // Number of indices

    // H5DWrite doesn't accept vectors so convert to array
    // **************************************************************    
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
    /*
        For writing k-wave float array to HDF5 file

        Input:      datasetName -- Name of variable
                    value       -- Array of values
    */

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
    /*
        Writes HDF5 file for k-Wave input for C++ binary for 3D data
    */

    std::cout << "Beginning to write K-Wave Input File." << std::endl;
        
    if(binDimSize != 3)
    {
        std::cout << "KWaveInput::WriteFile() only accepts 3D Cartesian, instead received a size: " << binDimSize << std::endl;
        return;
    }

    if(grueneisen != 1)
    {
        for(int i=0; i<binDim[0]*binDim[1]*binDim[2]; i++)
        {
            energyDensityData[i] *= grueneisen; 
        }
    }
    
    if( enableSmooth )
        fenergyDensityData = Smooth(energyDensityData, binDim, restoreMagnitude, rotWindow);
    else
    {
        fenergyDensityData = new float[binDim[0]*binDim[1]*binDim[2]]{0};
        std::copy(energyDensityData, energyDensityData+(binDim[0]*binDim[1]*binDim[2]), fenergyDensityData);
    }

    float gridSpacing[binDimSize];
    gridSpacing[0] = dx;
    gridSpacing[1] = dy;
    gridSpacing[2] = dz;

    Write1DFloat("dx",dx);
    Write1DFloat("dy",dy);
    Write1DFloat("dz",dz);

    float dt = CalculateDT(cfl, gridSpacing, mediumSoundSpeed);
    Write1DFloat("dt",dt);

    unsigned long Nt = CalculateNT(tEnd, dt);

    // Writing Variables to Input File
    // **************************************************************    
    Write1DLong("Nt",Nt);
    Write1DLong("Nx",Nx);
    Write1DLong("Ny",Ny);
    Write1DLong("Nz",Nz);

    Write1DLong("absorbing_flag",0);
    Write1DLong("axisymmetric_flag",0);
    
    Write1DFloat("c0",mediumSoundSpeed);
    Write1DFloat("c_ref",cRef);

    Write1DLong("elastic_flag",0);
    Write1DLong("nonlinear_flag",0);
    Write1DLong("nonuniform_grid_flag",0);

    Write1DLong("p0_source_flag",1);
    WritePData("p0_source_input", Nz, Ny, Nx, fenergyDensityData);
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

void KWaveInput::WriteFileAS()
{
    /*
        Writes HDF5 file for k-Wave input for C++ binary for 2D data (axisymmetric)
    */

    std::cout << "Beginning to write K-Wave Input File." << std::endl;

    if(binDimSize != 2)
    {
        std::cout << "KWaveInput::WriteFileAS() only accepts 2D Cartesian, instead received a size: " << binDimSize << std::endl;
        return;
    }
    
    if(grueneisen != 1)
    {
        for(int i=0; i<binDim[0]*binDim[1]; i++)
        {
            energyDensityData[i] *= grueneisen; 
        }
    }
    
    if( enableSmooth )
        fenergyDensityData = SmoothAS(energyDensityData, binDim, restoreMagnitude, rotWindow);
    else
    {        
        fenergyDensityData = new float[binDim[0]*binDim[1]]{0};
        std::copy(energyDensityData, energyDensityData+(binDim[0]*binDim[1]), fenergyDensityData);
    }

    float gridSpacing[binDimSize];
    gridSpacing[0] = dx;
    gridSpacing[1] = dy;

    Write1DFloat("dx",dx);
    Write1DFloat("dy",dy);

    float dt = CalculateDT(cfl, gridSpacing, mediumSoundSpeed);
    Write1DFloat("dt",dt);

    unsigned long Nt = CalculateNT(tEnd, dt);

    // Writing Variables to Input File
    // **************************************************************    
    Write1DLong("Nt",Nt);
    Write1DLong("Nx",Nx);
    Write1DLong("Ny",Ny);
    Write1DLong("Nz",Nz);

    Write1DLong("absorbing_flag",0);
    Write1DLong("axisymmetric_flag",1);
    
    Write1DFloat("c0",mediumSoundSpeed);
    Write1DFloat("c_ref",cRef);

    Write1DLong("elastic_flag",0);
    Write1DLong("nonlinear_flag",0);
    Write1DLong("nonuniform_grid_flag",0);

    
    Write1DLong("p0_source_flag",1);
    WritePData("p0_source_input", Nz, Ny, Nx, fenergyDensityData);
    Write1DLong("p_source_flag",0);
    
    
    Write1DFloat("pml_x_alpha",pmlAlphaX);
    Write1DLong("pml_x_size",pmlSizeX);
    Write1DFloat("pml_y_alpha",pmlAlphaY);
    Write1DLong("pml_y_size",pmlSizeY);
    Write1DLong("pml_z_size",pmlSizeZ);
    
    Write1DFloat("rho0",1);
    Write1DFloat("rho0_sgx",1);
    Write1DFloat("rho0_sgy",1);

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
