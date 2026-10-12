#include <Rcpp.h>
#include "AnisotropicDiffusion.h"
using namespace Rcpp;

// [[Rcpp::export]]
double fsd(SEXP x)
{
  double M = 0.0;
  double S = 0.0;
  int k = 1;
  int n = LENGTH(x);

  if (n == 1 )
    return NA_REAL;

  if (Rf_isReal(x))
  {
    double* px = REAL(x);

    for (int i = 0; i < n; ++i)
    {
      double tmpM = M;
      M += (px[i] - tmpM) / k;
      S += (px[i] - tmpM) * (px[i] - M);
      k++;
    }
  }
  else if (Rf_isInteger(x))
  {
    int* px = INTEGER(x);

    for (int i = 0; i < n; ++i)
    {
      double tmpM = M;
      M += (px[i] - tmpM) / k;
      S += (px[i] - tmpM) * (px[i] - M);
      k++;
    }
  }
  else
  {
    Rf_error("x is not INTEGER or REAL");
  }

  return std::sqrt(S / (k-2));
}

// [[Rcpp::export]]
IntegerMatrix anisotropic_diffusion(Rcpp::IntegerMatrix image, unsigned int iterations, double lambda, double k)
{
  AnisotropicDiffusion AD(Rcpp::clone(image), iterations, lambda, k);
  AD.applyDiffusion();
  return AD.getImage();
}
