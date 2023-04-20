#include "fonction.h"
#include <cmath>
#include <vector>

//////////////////////////////////////////////////////////////////////////
//// fonction f,g et h en fonction des différents cas
//////////////////////////////////////////////////////////////////////////

double f(double x, double y, double Lx, double Ly, int cas)
{
    // fonction f pour chacun des cas
    double res;
    if(cas==1)
    {
        res=2*(y-y*y+x-x*x);
    }
    if(cas==2)
    {
        res=sin(x)+cos(y);
    }
    if(cas==3)
    {
        res=exp(-pow(x-(Lx/2),2))+exp(-pow(y-(Ly/2),2));
    }
    return res;
}

double g(double x, double y, double Lx, double Ly, int cas)
{
    // fonction g pour chacun des cas
    if(cas==1)
    {
        return 0.0 ;
    }
    if(cas==2)
    {
        return sin(x)+cos(y);
    }
    if(cas==3)
    {
        return 0.0 ;
    }
}

double h(double x, double y, double Lx, double Ly, int cas)
{
    // fonction h pour chacun des cas
    if(cas==1)
    {
        return 0.0 ;
    }
    if(cas==2)
    {
        return sin(x)+cos(y);
    }
    if(cas==3)
    {
        return 1.0 ;
    }
}

vector<double> F_b(double Lx,double Ly, double alpha, double beta, double gamma, int cas, int Nx, int Ny, vector<double> U)
{
    // fonction pour avoir le vecteur F du second membre 
    vector <double> F;

    //F[0]=; // cas j=0 (j boucle dans les blocs)

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx-1
    {
        //F[j]=;
    }

      // cas j=Nx-1

    /////////////////////////////////////////////////////////////////
    /// cas pour i différent de 0 et Ny-1
    /////////////////////////////////////////////////////////////////
    for (int i=1;i<Ny-1;i++) 
    {
         // cas j=0

        for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
        {
            
        }

         // cas j=Nx-1
    }

    /////////////////////////////////////////////////////////////////
    /// cas i=Ny-1
    /////////////////////////////////////////////////////////////////
    

     // cas j=0

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
    {
        
    }

     // cas j=Nx-1

    return F;
}


//////////////////////////////////////////////////////////////////////////
//// fonction multiplication matrice A* vecteur U
//////////////////////////////////////////////////////////////////////////

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

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx-1
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
