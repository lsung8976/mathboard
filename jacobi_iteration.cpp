#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;


int main()
{
    std::ofstream out("jacobiTest.txt", std::ios::app);
    
    int max_iter = 200;
    double x[max_iter+2], y[max_iter+2], z[max_iter+2];

    x[0] = y[0] = z[0] = 0;

    for(int i = 1;i <= max_iter; i++)
    {
        x[i] = (3 - 2 * y[i-1] - z[i-1])/4; 
        y[i] = (-1 - x[i-1] - z[i-1])/3;
        z[i] = (4 - x[i-1] - y[i-1])/4;
        if (out.is_open()) {
                out << x[i] << " " << y[i] << " " << z[i] << endl;
        }
    }

    return 0;
}