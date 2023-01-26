#ifndef POLYFIT_H
#define POLYFIT_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <iomanip>

class PolyFit
{
public:
    PolyFit();
    ~PolyFit();
    double incbeta(double a, double b, double x);
    double invincbeta(double y, double alpha, double beta);
    double CalculateTValueStudent(const double nu, const double alpha);
    double cdfStudent(const double nu, const double t);
    double cdfFisher(const double df1, const double df2, const double x);
    double **Make2DArray(const size_t rows, const size_t cols);
    double **MatTrans(double **array, const size_t rows, const size_t cols);
    double **MatMul(const size_t m1, const size_t m2, const size_t m3, double **A, double **B);
    void MatVectMul(const size_t m1, const size_t m2, double **A, double *v, double *Av);
    double determinant(double **a, const size_t k);
    void transpose(double **num, double **fac, double **inverse, const size_t r);
    void cofactor(double **num, double **inverse, const size_t f);
    void displayMat(double **A, const size_t n, const size_t m);
    double CalculateRSS(const double *x, const double *y, const double *a, double **Weights,
                        const bool fixed, const size_t N, const size_t n);
    double CalculateTSS(const double *x, const double *y, const double *a, double **Weights, 
                        const bool fixed, const size_t N, const size_t n);
    double CalculateR2COD(const double *x, const double *y, const double *a, double **Weights,
                            const bool fixed, const size_t N, const size_t n);
    double CalculateR2Adj(const double *x, const double *y, const double *a, double **Weights,
                            const bool fixed,const size_t N, const size_t n);
    void Fit(const double *x, double *y, const size_t n, const size_t k, const bool fixedinter,
                const double fixedinterval, double *beta, double **Weights, double **XTWXInv);
    double calculatePoly(const double x, const double *a, const size_t n);
    void WriteCIBands(std::string filename, const double *x, const double *coefbeta, double **XTXInv, 
                        const double tstudentval, const double SE, const size_t n, const size_t k);
    void CalculateWeights(const double *erry, double **Weights, const size_t n,
                            const int type);
    void CalculateSERRBeta(const bool fixedinter, const double SE, size_t k, double *serbeta, double **XTWXInv);
    void DisplayPolynomial(const size_t k);
    void DisplayANOVA(const size_t nstar, const size_t k, const double TSS, const double RSS);
    void DisplayCoefs(const size_t k, const size_t nstar, const double tstudentval, const double *coefbeta, const double *serbeta);
    void DisplayStatistics(const size_t n, const size_t nstar, const size_t k, const double RSS, const double R2,
                            const double R2Adj, const double SE);
    void DisplayCovCorrMatrix(const size_t k, const double sigma, const bool fixed, double **XTWXInv);
    
};
#endif // POLYFIT_H
