#include "MatFile.h"

using namespace std;

#include <cstdio>
#include <iostream>
#include <cstring>

MatFile::MatFile() : file(-1)
{
}
MatFile::MatFile(const string& fn) : file(-1)
{
	Open(fn);
}

MatFile::~MatFile()
{
	Close();
}

bool MatFile::Open(const string& fn)
{
	if (file >= 0)
		Close();

	filename = fn;

	file = H5Fcreate(fn.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	return file >= 0;
}

bool MatFile::Close()
{
	if (file >= 0)
	{
		H5Fclose(file);
		file = -1;
		return PrependHeader(filename.c_str());
	}
	return false;
}

bool MatFile::IsOpen() const
{
	return file >= 0;
}

bool MatFile::Write3Dim(const string& varname, int nx, int ny, int nz, const double* data)
{
	if(file < 0)
		return false;

	if(nx <= 0 || ny <= 0 || nz <= 0)
		return true;

	hsize_t dimsf[3];
	dimsf[0] = nx;
	dimsf[1] = ny;
    dimsf[2] = nz;

	return H5LTmake_dataset_double(file, varname.c_str(), 3, dimsf, data) >= 0;
}

bool MatFile::WriteMatrix(const string& varname, int nx, int ny, const double* data)
{
	if (file < 0)
		return false;

	if (nx <= 0 || ny <= 0)
		return true;

	hsize_t dimsf[2];
	dimsf[0] = nx;
	dimsf[1] = ny;
	return H5LTmake_dataset_double(file, varname.c_str(), 2, dimsf, data) >= 0;
}

bool MatFile::WriteVector(const string& varname, int nx, const double* data)
{
	if (file < 0)
		return false;

	if (nx <= 0)
		return true;

	hsize_t dimsf[1];
	dimsf[0] = nx;
	return H5LTmake_dataset_double(file, varname.c_str(), 1, dimsf, data) >= 0;
}



bool MatFile::PrependHeader(const char* filename)
{
	// Prepend the header.
	FILE* fd = fopen(filename, "r+b");
	if (!fd)
	{
		cerr << "Couldn't open file for writing header." << endl;
		return false;
	}

	char header[512];
	memset(header, 0, sizeof(header));
	sprintf(header, "MATLAB 7.3 format file written by MatFile class, by Tim Hutt"); // Do not make this longer than 116 chars (I think)
	header[124] = 0;
	header[125] = 2;
	header[126] = 'I';
	header[127] = 'M';

	// Get file length.
	fseek(fd, 0, SEEK_END);
	long long length = ftell(fd);
	fseek(fd, 0, SEEK_SET);

	// TODO: Do this properly without reading entire file into memory.
	if (length > 1024L*1024L*1024L*10L) // 10 GB.
	{
		cerr << "File too big to write header, sorry.";
		fclose(fd);
		return false;
	}

	unsigned char* buffer = new unsigned char[length];
	if (fread(buffer, 1, length, fd) != length)
	{
		delete[] buffer;
		cerr << "Couldn't read file, sorry.";
		fclose(fd);
		return false;
	}

	fseek(fd, 0, SEEK_SET);

	fwrite(header, 1, sizeof(header), fd);
	fwrite(buffer, 1, length, fd);
	std::cout << " MatFile Written." << std::endl;
	fclose(fd);

	delete[] buffer;

	return true;
}
