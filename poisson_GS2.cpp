#include <iostream>
#include <math.h>
#include <fstream>
#define _USE_MATH_DEFINES
using namespace std;


int main()
{
    std::ofstream out("test22.txt", std::ios::app);
    std::ofstream exout("test23.txt", std::ios::app);
    
    int max_iter = 10000;
    int Nx = 10;
    int Ny = Nx;

    double x[Nx+2], y[Ny+2], f[Nx+2][Ny+2], u[Nx+2][Ny+2], exu[Nx+2][Ny+2];
    double src;

    double h = 1.0 / Nx;
    double h2 = h * h;

    double x_size =  (1.0 - h)/Nx;
    double y_size =  (1.0 - h)/Ny;

    double coef = 0;

    for(int i = 0 ;i <= Nx; i++)
    {
        x[i] = i * x_size;
        y[i] = i * y_size;
    }
    
    for(int i = 0 ;i <= Nx+1; i++)
    {
        for(int j = 0 ;j <= Ny+1; j++)
        {
            u[i][j] = 0.0;
            exu[i][j] = 0.5 * cos(2 * M_PI * x[i]) * cos(2 * M_PI * y[j]);
            f[i][j] = (-8) * M_PI * M_PI * exu[i][j];
        }
    }

    for(int k = 1; k <= max_iter; k++)
    {
        for(int i = 1; i <= Nx; i++)
        {
            for(int j = 1; j <= Ny; j++)
            {
                if(i != 1 || j != 1)
                {
                    src = - h2 * f[i][j]; coef = 0.0;
                    if(i>1){
                        src = src + u[i-1][j]; coef++;
                    }
                    if(i<Nx){
                        src = src + u[i+1][j]; coef++;
                    }
                    if(j>1){
                        src = src + u[i][j-1]; coef++;
                    }
                    if(j<Ny){
                        src = src + u[i][j+1]; coef++;
                    }
                
                    u[i][j] = src / coef;
                }
            }
        }
    }

    for(int i = 1; i <= Nx; i++)
    {
        for(int j = 1; j <= Ny; j++)
        {
            if (out.is_open() && exout.is_open())
            {
                out << x[i] << " " << y[j] << " " << u[i][j] << endl;
                exout << x[i] << " " << y[j] << " " << exu[i][j] << endl;
            }
        }
    }

    return 0;
}