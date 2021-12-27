#include <iostream>
#include <math.h>
#include <fstream>
#define _USE_MATH_DEFINES
using namespace std;


int main()
{
    std::ofstream out("poissonTest.txt", std::ios::app);
    
    int max_iter = 10;
    int N = 10;
    double x[N+2], f[N+2], t[10*N+2], u[N+2];
    
    double h = 1.0 / N;

    double x_size =  (1 - h)/N;
    double t_size = (1 - h)/ (10*N);
    for(int i = 1 ;i <= N; i++)
    {
        x[i] = i * x_size;
        t[i] = i * t_size;
        f[i] = -4 * M_PI * M_PI * cos(2*M_PI*x[i]);
    }
    for(int i = 1 ;i <= N; i++)
        u[i] = 0;
    
    for(int k = 1;k <= max_iter; k++)
    {
        for(int i = 1; i <= N-1; i++)
        {
            if(i>1)
                u[i] = (u[i+1]+u[i-1]-h*h*f[i])/2;
            if (out.is_open())
                out << k << " " << x[i] << " " << u[i] << endl;
        }
        u[N] = u[N-1]-h*h*f[N];
    }

    return 0;
}