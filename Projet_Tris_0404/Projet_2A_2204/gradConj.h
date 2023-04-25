
#include <vector>
using namespace std;
typedef vector<double> Vector;
double dot(const Vector& x, const Vector& y);
double norm(const Vector& x) ;
vector<double> gradientConj (Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b);