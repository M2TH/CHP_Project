


#include "GradConjMPI.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "fonction.h"
#include <mpi.h>

using namespace std;

typedef vector<double> Vector;

// Le produit scalaire des vecteurs x et y
double dot(const Vector& x, const Vector& y) {
    int n = x.size();
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

// La norme 2 d'un vecteur x
double norm(const Vector& x) {
    int n = x.size();
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += x[i] * x[i];
    }
    return sqrt(result);
}


// // Compute the dot product of two vectors on different processors and return the result
// double dot_mpi(const Vector& x, const Vector& y, int iBeg, int iEnd, int np, int me) {
//     int n = x.size();
//     double local_result = dot(x, y, iBeg, iEnd);

//     double result;
//     MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//     return result;
// }

// // Compute the norm of a vector on different processors and return the result
// double norm_mpi(const Vector& x, int iBeg, int iEnd, int np, int me) {
//     int n = x.size();
//     double local_result = norm(x, iBeg, iEnd);

//     double result;
//     MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

//     return result;
// }

// Apply the matrix A to a vector x on different processors and return the result
Vector matvec_mpi(double alpha, double beta, double gamma, int Nx, int Ny, const Vector& x, int iBeg, int iEnd, int np, int me) {
    int n = Nx * Ny;
    Vector AK(n);
    for (int i = iBeg; i <= iEnd; i++) {
        int ix = i % Nx;
        int iy = i / Nx;
        AK[i] = alpha * x[i] + beta * (ix == 0 ? 0 : x[i-1]) + beta * (ix == Nx-1 ? 0 : x[i+1]) + gamma * (iy == 0 ? 0 : x[i-Nx]) + gamma * (iy == Ny-1 ? 0 : x[i+Nx]);
    }

    Vector W(n);
    MPI_Allreduce(&AK[0], &W[0], n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); //dans AK[me] ? 

    return W;
}





// Solveur du system AU=F /Ax=b avec le gradient conjugué 
vector<double> gradientConj_mpi(Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b, const int iBeg, const int iEnd, const int np, const int me) { //(Vector& K, double tol, int maxiter, double alpha, double beta, double gamma, Vector& b) {    
    // tol est la tolérance (seuil à atteindre pour considere la convergence)
    // maxiter est le garde-fou pour arreter le calcul en cas de non convergence
    // K vecteur condition initiale     
    // double alpha, beta, gamma; // Inclu alpha, beta, gamma (les coefficient de la matrice) pris en paramètre 
    // b : vecteur second membre

    Vector AK = matvec_mpi(alpha, beta, gamma, Nx, Ny, K, iBeg, iEnd, np, me); // matvec(alpha, beta, gamma, Nx, Ny, K);
    
    int n = Nx* Ny ;

    Vector r(n), d(n) ;

    for (int i = 0; i < n; i++) {
        r[i] = AK[i] - b[i];
        d[i] = r[i];
    }
    
    // Initialisation
    Vector x = K;
    Vector W(n);
    
    //produit W = A*d
    // for (int i = 0; i < n; i++) {
    //     W[i]=0.0;
    // }
    W = matvec_mpi(alpha, beta, gamma, Nx, Ny, K, iBeg, iEnd, np, me);
    //W = matvec(alpha, beta, gamma, Nx, Ny, d);
    
    // Iteration jusqu'à convergence ou max d'itération atteind
    int iter = 0;
    double coef_a, coef_b, residual;
    residual = norm(r);
    
    while ( (iter < maxiter) && (residual > tol)) {
        //produit W = A*d
        //W = matvec(alpha, beta, gamma, Nx, Ny, d);
        W = matvec_mpi(alpha, beta, gamma, Nx, Ny, K, iBeg, iEnd, np, me);

        // coef_a = (r*r) / (d*W)
        coef_a = dot(d, r) / dot(d, W);

        // Màj solution x et residu r
        for (int i = 0; i < n; i++) {
            x[i] -= coef_a * d[i];
            r[i] -= coef_a * W[i];
        }
        
        // Verif convergence
        residual = norm(r);
        if (residual < tol) {
            cout << "Convergence reached after " << iter << " iterations" << endl;
            break;
        }

        // Calcul coefficient beta du grad conj (pas la beta de la matrice)
        coef_b = dot(r, r) / dot(d, W);
        
        // Actualisation de la direction de recherhce
        for (int i = 0; i < n; i++) {
            d[i] = r[i] + coef_b * d[i];
        }
    
        iter++;
    }

    return x;
}





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ------------------------------ARCHIVES-------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// fonction charge qui permettait de tester le produit matvec sans utiliser une autre fonction

// void charge(int me, int n, int np, int *iBeg, int *iEnd)
// {
//     int divi, rest;

//     rest=n%np;
//     divi=n/np;
//     if(me<rest)
//     {
//         *iBeg=me*divi+me;
//         *iEnd=(me+1)*divi+me;
//     }
//     else
//     {
//         *iBeg=me*divi+rest;
//         *iEnd=(me+1)*divi-1+rest;
//     }
// }



// /// main servant à tester le produit matrice vecteur parallélisé
// int main(int argc, char **argv)
// {

//     int me,nproc,n;

//     MPI_Init(&argc,&argv);
//     MPI_Comm_size(MPI_COMM_WORLD,&nproc);
//     MPI_Comm_rank(MPI_COMM_WORLD,&me);

//     Vector X(n) ; 
//     for(int i=0;i<n;i++)  //vecteur x du produit Ax=b
//     {
//         X[i]=1.0 ;
//     }

//     double alpha, beta, gamma;
//     alpha=1. ;
//     beta = 1. ;
//     gamma = 1. ;

//     int Nx;
//     int Ny; 
//     Nx = 3; // matrice de taille n*n (test simple)
//     Ny = 3; 
//     n= Nx*Ny;

//     int iBeg,iEnd;
//     charge(me,n,nproc,&iBeg,&iEnd);
//     printf("Je suis le proc %d avec iBeg=%d et iEnd=%d\n",me,iBeg,iEnd);

//     Vector W(n);
//     W = matvec_mpi( alpha, beta, gamma, Nx, Ny, X , iBeg, iEnd, nproc, me);

//     if(me==0){
//         printf("Vecteur solution B \n");
//         for(int i=0;i<n;i++){
//             printf(" %f  \n", W[i]) ; 
//         }
//     }
    
//     MPI_Finalize() ; 

//     return 0; 

// }


