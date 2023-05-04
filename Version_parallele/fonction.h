#include <vector>
using namespace std;
double f(double x, double y, double t, double Lx, double Ly, int cas);
double g(double x, double y, double Lx, double Ly, int cas);
double h(double x, double y, double Lx, double Ly, int cas);
vector <double> F_b(double Lx,double Ly, double beta, double gamma, double dt, double t, double D, int cas, int Nx, int Ny, vector<double> U);
// vector <double> matvec(double alpha, double beta, double gamma, int Nx, int Ny, vector<double> U);
void charge(int me, int n, int np, int *iBeg, int *iEnd);
