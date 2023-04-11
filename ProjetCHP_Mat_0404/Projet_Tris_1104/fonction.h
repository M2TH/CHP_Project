#include <vector>
using namespace std;
double f_cas1(double x, double y, double t, double Lx, double Ly);
double g_cas1(double x, double y, double t, double Lx, double Ly);
double h_cas1(double x, double y, double t, double Lx, double Ly);
double f_cas2(double x, double y, double t, double Lx, double Ly);
double g_cas2(double x, double y, double t, double Lx, double Ly);
double h_cas2(double x, double y, double t, double Lx, double Ly);
double f_cas3(double x, double y, double t, double Lx, double Ly);
double g_cas3(double x, double y, double t, double Lx, double Ly);
double h_cas3(double x, double y, double t, double Lx, double Ly);
vector <double> matvec(double alpha, double beta, double gamma, int Nx, int Ny, vector<double> U);
void charge(int me, int n, int np, int *iBeg, int *iEnd);