#include "function.h"
#include <cmath>
#include <vector>

//////////////////////////////////////////////////////////////////////////
//// fonction cas 1
//////////////////////////////////////////////////////////////////////////

double f_cas1(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=2*(y-y*y+x-x*x);
    
    return res;
}

double g_cas1(double x, double y, double t, double Lx, double Ly)
{
    return 0;
}

double h_cas1(double x, double y, double t, double Lx, double Ly)
{
    return 0;
}

//////////////////////////////////////////////////////////////////////////
//// fonction cas 2
//////////////////////////////////////////////////////////////////////////

double f_cas2(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=sin(x)+cos(y);
    
    return res;
}

double g_cas2(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=sin(x)+cos(y);
    
    return res;
}

double h_cas2(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=sin(x)+cos(y);
    
    return res;
}

//////////////////////////////////////////////////////////////////////////
//// fonction cas 3
//////////////////////////////////////////////////////////////////////////

double f_cas3(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=exp(-pow(x-(Lx/2),2))+exp(-pow(y-(Ly/2),2));
    
    return res;
}

double g_cas3(double x, double y, double t, double Lx, double Ly)
{
    return 0;
}

double h_cas3(double x, double y, double t, double Lx, double Ly)
{
    return 1;
}

vector<double> matvec(double alpha, double beta, double gamma, int Nx, int Ny, vector<double> U)
{
    vector<double> res(Nx*Ny);
    // alpha : diagonale central
    // beta : surdiagonale
    // gamma : diagonale externe
    /////////////////////////////////////////////////////////////////
    /// cas i=0 (i boucle sur les blocs)
    /////////////////////////////////////////////////////////////////

    res[0]=alpha*U[0]+beta*U[1]+gamma*U[Nx]; // cas j=0 (j boucle dans les blocs)

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
    {
        res[j]=beta*U[j-1]+alpha*U[j]+beta*U[j+1]+gamma*U[j+Nx];
    }

    res[Nx-1]=beta*U[Nx-2]+alpha*U[Nx-1]+gamma*U[Nx+Nx-1]; // cas j=Nx-1

    /////////////////////////////////////////////////////////////////
    /// cas pour i différent de 0 et Ny-1
    /////////////////////////////////////////////////////////////////
    for (int i=1;i<Ny-1;i++) 
    {
        int dy; // indice du alpha en haut à gauche du bloc
        dy=i*Nx;
        res[dy]=gamma*U[dy-Nx]+alpha*U[dy]+beta*U[dy+1]+gamma*U[dy+Nx]; // cas j=0

        for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
        {
            res[j+dy]=gamma*U[dy-Nx+j]+beta*U[dy+j-1]+alpha*U[dy+j]+beta*U[dy+j+1]+gamma*U[dy+j+Nx];
        }

        res[dy+Nx-1]=gamma*U[dy-1]+beta*U[dy+Nx-2]+alpha*U[dy+Nx-1]+gamma*U[dy+Nx]; // cas j=Nx-1
    }

    /////////////////////////////////////////////////////////////////
    /// cas i=Ny-1
    /////////////////////////////////////////////////////////////////
    int dy;
    dy=Nx*(Ny-1);

    res[dy]=gamma*U[dy-Nx]+alpha*U[dy]+beta*U[dy+1]; // cas j=0

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
    {
        res[dy+j]=gamma*U[dy+j-Nx]+beta*U[dy+j-1]+alpha*U[dy+j]+beta*U[dy+j+1];
    }

    res[dy+Nx-1]=gamma*U[dy-1]+beta*U[dy+Nx-2]+alpha*U[dy+Nx-1]; // cas j=Nx-1

    return res;
}

void charge(int me, int n, int np, int *iBeg, int *iEnd)
{
    int divi, rest;

    rest=n%np;
    divi=n/np;
    if(me<rest)
    {
        *iBeg=me*divi+me;
        *iEnd=(me+1)*divi+me;
    }
    else
    {
        *iBeg=me*divi+rest;
        *iEnd=(me+1)*divi-1+rest;
    }
}
