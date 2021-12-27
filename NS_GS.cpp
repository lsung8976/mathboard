#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>

using namespace std;

int main()
{  
    std::ofstream out("test66.txt", std::ios::app);

    int nx = 16;
    int ny = nx;
    int max_it = 10000000;
    
    double xleft = 0.0;
    double xright = 1.0;
    double yleft = 0.0;
    double yright = 1.0;

    double h = 1 / double(nx);
    double h2 = h * h;

    double dt = 0.001;
    double Re = 100.0;
    double a = 1.0;
    double tol, coef, src;
    double resid;

    double x_size, y_size;

    double u[nx+2][ny+3], nu[nx+2][ny+3], v[nx+3][ny+2], nv[nx+3][ny+2], p[nx+3][ny+3], oldp[nx+3][ny+3],  adv_u[nx+3][ny+3], adv_v[nx+3][ny+3]
                , tu[nx+2][ny+2], tv[nx+2][ny+2], x[nx+2], y[ny+2], vel_u[nx+2][ny+3], vel_v[nx+3][ny+2];
    double dpdx[nx+2][ny+2], dpdy[nx+2][ny+2];

    for(int i = 1; i <= nx+2; i++)
    {
        for(int j = 1; j<=ny+2; j++)
        {
            if(j < ny+2)
            {
                v[i][j] = 0.0;
                nv[i][j] = 0.0;
            }
            if(i < nx+2)
            {
                u[i][j] = 0.0;
                nu[i][j] = 0.0;
            }
            p[i][j] = 0.0;
        }
    }
    
    
    for(int it = 1; it <= max_it; it++)
    {
        //it 돌때마다 조건 초기화
        for(int j = 2;j <= ny+1 ;j++)
        {
            u[1][j]= u[nx+1][j] = 0;
        }
        for(int i = 2;i <= nx+1 ;i++)
        {
            v[i][1] = v[i][ny+1] = 0.0;
        }
        for(int i = 1;i <= nx+1; i++)
        {
            u[i][1] = -u[i][2];
            u[i][ny+2] = 2 * a -u[i][ny+1];
        }
        for(int j = 1;j <= ny+1; j++)
        {
            v[1][j] = -v[2][j];
            v[nx+2][j] = -v[nx+1][j];
        }

        for(int i = 2;i <= nx+1 ;i++)
        {
            adv_u[1][i] = adv_u[nx+1][i] = 0.0;
        }

        for(int i=2;i<=nx;i++)
        {
            for(int j=2;j<=ny+1;j++)
            {
                if(j==2)
                {
                    adv_u[i][j] = 0.5 * u[i][j] * (u[i+1][j]-u[i-1][j])/h
                                    + 0.25 * (v[i][j] + v[i+1][j]) * 0.5 * (u[i][j+1]-(-u[i][j]))/h;
                }
                else if(i==ny+1)
                {
                    adv_u[i][j] = 0.5 * u[i][j] * (u[i+1][j]-u[i-1][j])/h
                                    + 0.25 * (v[i][j-1] + v[i+1][j-1]) * 0.5 * ((2 * a -u[i][j])-(u[i][j-1]))/h;
                }
                else
                {
                    adv_u[i][j] = 0.5 * u[i][j] * (u[i+1][j]-u[i-1][j])/h
                                    + 0.25 * (v[i][j-1] + v[i+1][j-1]) + (v[i][j] + v[i+1][j]) * 0.5 * (u[i][j+1] - u[i][j-1])/h;
                }
            }
        }

        for(int i = 2;i <= nx+1 ;i++)
        {
            adv_v[i][1] = adv_v[i][ny+1] = 0.0;
        }

        for(int i=2;i<=nx + 1;i++)
        {
            for(int j=2;j<=ny;j++)
            {
                if(j==2)
                {
                    adv_v[i][j] = 0.5 * v[i][j] * (v[i][j+1]-v[i][j-1])/h
                                    + 0.25 * (u[i][j] + u[i][j+1]) * 0.5 * (v[i+1][j]-(-v[i][j]))/h;
                }
                else if(i==ny+1)
                {
                    adv_v[i][j] = 0.5 * v[i][j] * (v[i][j+1]-v[i][j-1])/h
                                    + 0.25 * (u[i-1][j] + u[i-1][j+1]) * 0.5 * (-v[i][j]-v[i-1][j])/h;
                }
                else
                {
                    adv_v[i][j] = 0.5 * v[i][j] * (v[i][j+1]-v[i][j-1])/h
                                    + 0.25 * (u[i-1][j] + u[i][j] + u[i-1][j+1] + u[i][j+1])
                                        * 0.5 * (v[i+1][j] - v[i-1][j])/h;
                }
            }
        }

        for(int i=2;i<=nx;i++)
        {
            for(int j=2;j<=ny+1;j++)
            {
                tu[i][j] = u[i][j] + dt * ((u[i+1][j] + u[i-1][j] -4.0*u[i][j] + u[i][j+1] + u[i][j-1]) / (Re * h2) - adv_u[i][j] );
            }
        }

        for(int j=2;j<=ny+1;j++)
            tu[1][j] = tu[nx+1][j] = 0.0;

        for(int i=2;i<=nx+1;i++)
        {
            for(int j=2;j<=ny;j++)
            {
                tv[i][j] = v[i][j] + dt * ((v[i+1][j] + v[i-1][j] -4.0*v[i][j] + v[i][j+1] + v[i][j-1]) / (Re * h2) - adv_v[i][j] );
            }
        }
        
        for(int i=2; i<=nx+1; i++)
            tv[i][1] = tv[i][ny+1] = 0.0;

        tol = 1.0e-4;
        resid = 1.0;
        
        for(int i=1;i<=nx+2;i++)
        {
            for(int j=1;j<=ny+2;j++)
            {
                oldp[i][j] = p[i][j];
            }
        }

        while (resid > tol)
        {
            for(int i=1;i<=nx+2;i++)
            {
                for(int j=1;j<=ny+2;j++)
                {
                    oldp[i][j] = p[i][j];
                }
            }

            for(int i=2;i<=nx+1;i++)
            {
                for(int j=2;j<=ny+1;j++)
                {
                    coef = 0.0;
                    src = - h2 * (tu[i][j]-tu[i-1][j]+tv[i][j]-tv[i][j-1])/(h*dt);
                    if(i>2)
                    {
                        coef = coef + 1.0;
                        src = src + p[i-1][j];
                    }
                    if(i<nx+1)
                    {
                        coef = coef + 1.0;
                        src = src + p[i+1][j];
                    }
                    if(j>2)
                    {
                        coef = coef + 1.0;
                        src = src + p[i][j-1]; 
                    }
                    if(j<ny+1)
                    {
                        coef = coef + 1.0;
                        src = src + p[i][j+1];
                    }
                
                    p[i][j] = src / coef;
                }
            }

            for(int i=1;i<=nx+2;i++)
            {
                for(int j=1;j<ny+2;j++)
                {
                    if(i==1 && j==1)
                        resid = abs(p[i][j]-oldp[i][j]);
                    resid = max(abs(p[i][j]-oldp[i][j]), resid);
                }
            }

        }

        for(int i=1;i<=nx+1;i++)
        {
            for(int j=1;j<=ny;j++)
            {
                if(i==1)
                    dpdx[i][j] = 0.0;
                else if(i==nx+1)
                    dpdx[i][j] = 0.0;
                else
                    dpdx[i][j] = (p[i+1][j+1] - p[i][j+1]) / h;
            }
        }

        for(int i=1;i<=nx;i++)
        {
            for(int j=1;j<ny+1;j++)
            {
                if(j==1)
                    dpdy[i][j] = 0.0;
                else if(j==ny+1)
                    dpdy[i][j] = 0.0;
                else
                    dpdy[i][j] = (p[i+1][j+1] - p[i+1][j]) / h;
            }
        }

        for(int i = 2;i <= nx + 1; i++)
        {
            for(int j = 2; j <= ny + 1;j++)
            {
                nu[i][j] = tu[i][j] - dt*dpdx[i][j-1];
                nv[i][j] = tv[i][j] - dt*dpdy[i-1][j];
            }
        }

        for(int i=1;i <= nx+2; i++)
        {
            nv[i][1] = nv[i][ny+1] = 0.0;
        }
        
        for(int i=1;i <= nx+2; i++)
        {
            for(int j=1;j <=ny+2; j++)
            {
                if(j < ny+2)
                    v[i][j] = nv[i][j];
                if(i < nx+2)
                    u[i][j] = nu[i][j];
            }
            
        }

    }

    x_size = (1.0 - h) / nx;
    y_size = (1.0 - h) / ny;

    for(int i=1;i <= nx; i++)
    {
        x[i] = i * x_size;
        y[i] = i * y_size;
    }
    
    
    for(int i = 1; i <= nx; i++)
    {
        for(int j = 1; j <= ny; j++)
        {
            vel_u[i][j] = (nu[i][j+1]+nu[i+1][j+1])/2;
            vel_v[i][j] = (nv[i+1][j]+nu[i+1][j+1])/2;
        }
    }

    for(int i = 1; i <= nx; i++)
    {
        for(int j = 1; j <= ny; j++)
        {
            if (out.is_open())
            {
                out << x[i] << " " << y[j] << " " << vel_u[i][j] << " " << vel_v[i][j] << endl;
            }
        }
    }


    return 0;
}