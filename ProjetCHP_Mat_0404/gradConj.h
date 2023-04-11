#include <vector>
using namespace std;
double dot(const Vector& x, const Vector& y);
double norm(const Vector& x) ;
void gradientConj(const Matrix& A, const Vector& F, Vector& U, double tol, int maxiter);


