#include <vector>
using namespace std;
typedef vector<double> Vector;
double dot(const Vector& x, const Vector& y);
double norm(const Vector& x) ;
Vector matvec_mpi(double alpha, double beta, double gamma, int Nx, int Ny, const Vector& x, int iBeg, int iEnd, int np, int me) ;
vector<double> gradientConj_mpi(Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b, const int iBeg, const int iEnd, const int np, const int me);
// vector<double> gradientConj(Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b);

