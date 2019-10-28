#include <iostream>
#include <fstream>
#include <cmath>

using fptr = double(double);  
int bisection(double xl, double xu, double Delta, fptr f);
int regulafalsi(double xl, double xu, double Delta, fptr f);
int newton_raphson(double x0, double Delta, fptr f, fptr fderiv);
double fun(double x);
double gaux(double x);
double deriv(double x);

int main(void)
{
 std::ofstream fout("1.txt");
 fout.precision(15); fout.setf(std::ios::scientific);
 fout << "# EPS-Delta" << "\t ";
 fout << "NBisect" << "\t";
 fout << "NFalsePos" << "\t";
 fout << "NNewton" << "\n";
  for (int ii=1; ii<= 9; ++ii)
    {
  double Delta = pow(10,-ii);
  double N = 0;
  fout << Delta << "\t  ";
  fout << bisection(0.1, 2.0, Delta, fun) <<"\t   ";
  fout << regulafalsi(0.1, 2.0, Delta, fun) << "\t ";
  fout << newton_raphson(0.1, Delta, fun, deriv) << "\n";
   }
  fout.close();
  return 0;
}

int bisection(double xl, double xu, double Delta, fptr f)
{
  double xr = xl;
  int N = 0;
  while(1) {
    N++;
    xr=(xl+xu)/2;
    if (std::fabs(f(xr)) <= Delta) {
      break;
    }
    else if (f(xr)*f(xl) < 0) {
      xu = xr;
    }
    else {
      xl = xr;
    }
 
  }

  return N;
}


int regulafalsi(double xl, double xu, double Delta, fptr f)
{
  double xr = xl;
  int N = 0;
  while(1) {
    N++;
    xr = xu - (f(xu)*(xl-xu))/(f(xl) -f(xu));
    if (std::fabs(f(xr)) <= Delta) {
      break;
    }
    else if (f(xr)*f(xl) < 0) {
      xu = xr;
    }
    else {
      xl = xr;
    }
    
  }


  return N;
}



double fun(double x)
{
  return 3*(1/exp(x)) - x;
}

double gaux(double x)
{
  return 3*(1/exp(x)) - x;
}
int newton_raphson(double x0, double Delta, fptr f, fptr fderiv)
{
  double xr = x0;
  int N = 0;
  while(1) {
    if (std::fabs(f(xr)) <= Delta) {
      break;
    }
    xr = xr - f(xr)/fderiv(xr);
    N++;
  }
  
  return N;
}


double deriv(double x)
{
  return -3*(1/exp(x)) - 1;
}
