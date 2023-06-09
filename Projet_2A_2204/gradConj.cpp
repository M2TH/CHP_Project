#include "gradConj.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "fonction.h"

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


// Solveur du system AU=F /Ax=b avec le gradient conjugué 
vector<double> gradientConj(Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b) { //(Vector& K, double tol, int maxiter, double alpha, double beta, double gamma, Vector& b) {    
    // tol est la tolérance (seuil à atteindre pour considere la convergence)
    // maxiter est le garde-fou pour arreter le calcul en cas de non convergence
    // K vecteur condition initiale     
    // double alpha, beta, gamma; // Inclu alpha, beta, gamma (les coefficient de la matrice) pris en paramètre 
    // b : vecteur second membre
    

    Vector AK = matvec(alpha, beta, gamma, Nx, Ny, K);
    
    // int n = A.size();
    int n = Nx* Ny ;

    Vector r(n), d(n), rkm1(n) ;

    // Initialisatioin residu r = b - A*x
    
    for (int i = 0; i < n; i++) {
        r[i] = AK[i] - b[i];
        d[i] = r[i];
    }
    
    // Initialisation
    Vector x = K;
    Vector W(n);
    
    //produit W = A*d
    for (int i = 0; i < n; i++) {
        W[i]=0.0;
    }
    // W = matvec(alpha, beta, gamma, Nx, Ny, d);
    
    // Iteration jusqu'à convergence ou max d'itération atteind
    int iter = 0;
    double coef_a, coef_b, residual;
    residual = norm(r);
    
    while ( (iter < maxiter) && (residual > tol)) {
        //produit W = A*d
        W = matvec(alpha, beta, gamma, Nx, Ny, d);

        // coef_a = (r*r) / (d*W)
        coef_a = dot(r, r) / dot(d, W);

        // Màj solution x er residu r
        for (int i = 0; i < n; i++) {
            x[i] -= coef_a * d[i];
            rkm1[i] = r[i]; 
            r[i] -= coef_a * W[i];
        }
        
        // Calcul coefficient beta du grad conj (pas la beta de la matrice)
        coef_b = dot(r, r) / dot(rkm1, rkm1);

        // Actualisation de la direction de recherhce
        for (int i = 0; i < n; i++) {
            d[i] = r[i] + coef_b * d[i];
        }
    
        // Verif convergence
        residual = norm(r);
        if (residual < tol) {
            cout << "Convergence reached after " << iter << " iterations" << endl;
            break;
        }

        iter++;
    }

    return x;
}




// vector<double> gradientConj(Vector& K, double tol, int maxiter, int Nx, int Ny, double alpha, double beta, double gamma, Vector& b) { //(Vector& K, double tol, int maxiter, double alpha, double beta, double gamma, Vector& b) {    
//     // tol est la tolérance (seuil à atteindre pour considere la convergence)
//     // maxiter est le garde-fou pour arreter le calcul en cas de non convergence
//     // K vecteur condition initiale     
//     // double alpha, beta, gamma; // Inclu alpha, beta, gamma (les coefficient de la matrice) pris en paramètre 
//     // b : vecteur second membre
    

//     Vector AK = matvec(alpha, beta, gamma, Nx, Ny, K);
    
//     // int n = A.size();
//     int n = Nx* Ny ;

//     Vector r(n), d(n) ;

//     // Initialisatioin residu r = b - A*x
    
//     for (int i = 0; i < n; i++) {
//         r[i] = AK[i] - b[i];
//         d[i] = r[i];
//     }
    
//     // Initialisation
//     Vector x = K;
//     Vector W(n);
    
//     //produit W = A*d
//     for (int i = 0; i < n; i++) {
//         W[i]=0.0;
//     }
//     W = matvec(alpha, beta, gamma, Nx, Ny, d);
    
//     // Iteration jusqu'à convergence ou max d'itération atteind
//     int iter = 0;
//     double coef_a, coef_b, residual;
//     residual = norm(r);
    
//     while ( (iter < maxiter) && (residual > tol)) {
//         //produit W = A*d
//         W = matvec(alpha, beta, gamma, Nx, Ny, d);

//         // coef_a = (r*r) / (d*W)
//         coef_a = dot(d, r) / dot(d, W);

//         // Màj solution x er residu r
//         for (int i = 0; i < n; i++) {
//             x[i] -= coef_a * d[i];
//             r[i] -= coef_a * W[i];
//         }
        
//         // Verif convergence
//         residual = norm(r);
//         if (residual < tol) {
//             cout << "Convergence reached after " << iter << " iterations" << endl;
//             break;
//         }

//         // Calcul coefficient beta du grad conj (pas la beta de la matrice)
//         coef_b = dot(r, r) / dot(d, W);
        
//         // Actualisation de la direction de recherhce
//         for (int i = 0; i < n; i++) {
//             d[i] = r[i] + coef_b * d[i];
//         }
    
//         iter++;
//     }

//     return x;
// }