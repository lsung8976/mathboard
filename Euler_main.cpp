#include<iostream>

#define nx 26//26
#define ny 26//26
#define nz 26//26
#define gamma 1.667
#define iter 5000

using namespace std;

int main()
{

int i, j, k, n, it; //iter;

double rho[nx+2][ny+2][nz+2]={0.0}, p[nx+2][ny+2][nz+2]={0.0};
double ux[nx+2][ny+2][nz+2]={0.0}, uy[nx+2][ny+2][nz+2]={0.0}, uz[nx+2][ny+2][nz+2]={0.0};
double ax[nx+2][ny+2][nz+2]={0.0}, ay[nx+2][ny+2][nz+2]={0.0}, az[nx+2][ny+2][nz+2]={0.0};
double bx[nx+2][ny+2][nz+2]={0.0}, by[nx+2][ny+2][nz+2]={0.0}, bz[nx+2][ny+2][nz+2]={0.0};

double it_inc=0.1, it_max=500.0, time=0.0; //426;
//iter=it_max/it_inc;


double ux_net[iter+1]={0.0}, uy_net[iter+1]={0.0}, uz_net[iter+1]={0.0}, u_net[iter+1]={0.0};
double bx_net[iter+1]={0.0}, by_net[iter+1]={0.0}, bz_net[iter+1]={0.0}, b_net[iter+1]={0.0};


double uxtemp=0.0, uytemp=0.0,uztemp=0.0, bxtemp=0.0, bytemp=0.0, bztemp=0.0;

/*
double k_rho[nx+2][ny+2][nz+2]={0.0}, k_p[nx+2][ny+2][nz+2]={0.0};
double k_ux[nx+2][ny+2][nz+2]={0.0}, k_uy[nx+2][ny+2][nz+2]={0.0}, k_uz[nx+2][ny+2][nz+2]={0.0};
double k_ax[nx+2][ny+2][nz+2]={0.0}, k_ay[nx+2][ny+2][nz+2]={0.0}, k_az[nx+2][ny+2][nz+2]={0.0};

double k1_rho[nx+2][ny+2][nz+2]={0.0}, k2_rho[nx+2][ny+2][nz+2]={0.0}, k3_rho[nx+2][ny+2][nz+2]={0.0}, k4_rho[nx+2][ny+2][nz+2]={0.0};
double k1_p[nx+2][ny+2][nz+2]={0.0}, k2_p[nx+2][ny+2][nz+2]={0.0}, k3_p[nx+2][ny+2][nz+2]={0.0}, k4_p[nx+2][ny+2][nz+2]={0.0};

double k1_ux[nx+2][ny+2][nz+2]={0.0}, k2_ux[nx+2][ny+2][nz+2]={0.0}, k3_ux[nx+2][ny+2][nz+2]={0.0}, k4_ux[nx+2][ny+2][nz+2]={0.0};
double k1_uy[nx+2][ny+2][nz+2]={0.0}, k2_uy[nx+2][ny+2][nz+2]={0.0}, k3_uy[nx+2][ny+2][nz+2]={0.0}, k4_uy[nx+2][ny+2][nz+2]={0.0};
double k1_uz[nx+2][ny+2][nz+2]={0.0}, k2_uz[nx+2][ny+2][nz+2]={0.0}, k3_uz[nx+2][ny+2][nz+2]={0.0}, k4_uz[nx+2][ny+2][nz+2]={0.0};

double k1_ax[nx+2][ny+2][nz+2]={0.0}, k2_ax[nx+2][ny+2][nz+2]={0.0}, k3_ax[nx+2][ny+2][nz+2]={0.0}, k4_ax[nx+2][ny+2][nz+2]={0.0};
double k1_ay[nx+2][ny+2][nz+2]={0.0}, k2_ay[nx+2][ny+2][nz+2]={0.0}, k3_ay[nx+2][ny+2][nz+2]={0.0}, k4_ay[nx+2][ny+2][nz+2]={0.0};
double k1_az[nx+2][ny+2][nz+2]={0.0}, k2_az[nx+2][ny+2][nz+2]={0.0}, k3_az[nx+2][ny+2][nz+2]={0.0}, k4_az[nx+2][ny+2][nz+2]={0.0};
*/

double mu=1.0, nu=0.006, eta=0.002, dt=0.01, t=0.0;
double dx=1.0, dy=1.0, dz=1.0;
double fvx=0.001, fvy=0.0, fvz=0.0;


FILE *fp = NULL;


for(i=0; i<=nx+1; i++){
       for(j=0; j<=ny+1; j++){
              for(k=0; k<=nz+1; k++){
                     rho[i][j][k]=0.1;
                     p[i][j][k]=0.00;
                     ux[i][j][k]=0.0;
                     uy[i][j][k]=0.0;
                     uz[i][j][k]=0.0;
                     ax[i][j][k]=0.01;
                     ay[i][j][k]=0.015;
                     az[i][j][k]=0.02;
              }
       }
}
//cout<<"No  x     y[0]      y[1]"<<endl;

printf("  it,  time,  i,   j,   k,  uxtemp,  uytemp,  uztemp,  bxtemp,  bytemp,  bztemp\n");



for(it=0; it<=iter; ++it){
       for(i=1; i<=nx; ++i){
              for(j=1; j<=ny; ++j){

                     rho[i][j][nz+1]=rho[i][j][1];
                     rho[i][j][0]=rho[i][j][nz]; //???????

                     p[i][j][nz+1]=p[i][j][1];
                     p[i][j][0]=p[i][j][nz];

                     ux[i][j][nz+1]=ux[i][j][1];
                     ux[i][j][0]=ux[i][j][nz];

                     uy[i][j][nz+1]=uy[i][j][1];
                     uy[i][j][0]=uy[i][j][nz];

                     uz[i][j][nz+1]=uz[i][j][1];
                     uz[i][j][0]=uz[i][j][nz];

                     ax[i][j][nz+1]=ax[i][j][1];
                     ax[i][j][0]=ax[i][j][nz];

                     ay[i][j][nz+1]=ay[i][j][1];
                     ay[i][j][0]=ay[i][j][nz];

                     az[i][j][nz+1]=az[i][j][1];
                     az[i][j][0]=az[i][j][nz];

              }
       }

       for(i=1; i<=nx; ++i){
              for(k=1; k<=nz; ++k){

                     rho[i][ny+1][k]=rho[i][1][k];
                     rho[i][0][k]=rho[i][ny][k];

                     p[i][ny+1][k]=p[i][1][k];
                     p[i][0][k]=p[i][ny][k];

                     ux[i][ny+1][k]=ux[i][1][k];
                     ux[i][0][k]=ux[i][ny][k];

                     uy[i][ny+1][k]=uy[i][1][k];
                     uy[i][0][k]=uy[i][ny][k];

                     uz[i][ny+1][k]=uz[i][1][k];
                     uz[i][0][k]=uz[i][ny][k];

                     ax[i][ny+1][k]=ax[i][1][k];
                     ax[i][0][k]=ax[i][ny][k];

                     ay[i][ny+1][k]=ay[i][1][k];
                     ay[i][0][k]=ay[i][ny][k];

                     az[i][ny+1][k]=az[i][1][k];
                     az[i][0][k]=az[i][ny][k];

              }
       }



       for(j=1; j<=nx; ++j){
              for(k=1; k<=nz; ++k){
                     rho[nx+1][j][k]=rho[1][j][k];
                     rho[0][j][k]=rho[nx][j][k];

                     p[nx+1][j][k]=p[1][j][k];
                     p[0][j][k]=p[nx][j][k];

                     ux[nx+1][j][k]=ux[1][j][k];
                     ux[0][j][k]=ux[nx][j][k];

                     uy[nx+1][j][k]=uy[1][j][k];
                     uy[0][j][k]=uy[nx][j][k];

                     uz[nx+1][j][k]=uz[1][j][k];
                     uz[0][j][k]=uz[nx][j][k];

                     ax[nx+1][j][k]=ax[1][j][k];
                     ax[0][j][k]=ax[nx][j][k];

                     ay[nx+1][j][k]=ay[1][j][k];
                     ay[0][j][k]=ay[nx][j][k];

                     az[nx+1][j][k]=az[1][j][k];
                     az[0][j][k]=az[nx][j][k];

              }
       }



       for(i=1; i<=nx; ++i){
              for(j=1; j<=ny; ++j){
                     for(k=1; k<=nz; ++k){
                            bx[i][j][k]=0.5*(az[i][j+1][k]-az[i][j-1][k])/dy - 0.5*(ay[i][j][k+1]-ay[i][j][k-1])/dz;
                            by[i][j][k]=0.5*(ax[i][j][k+1]-ax[i][j][k-1])/dz - 0.5*(az[i+1][j][k]-az[i-1][j][k])/dx;
                            bz[i][j][k]=0.5*(ay[i+1][j][k]-ay[i-1][j][k])/dx - 0.5*(ax[i][j+1][k]-ax[i][j-1][k])/dy;
                            ////////cout<<"it, i, j, k, bx, by, bz = "<<it<<' '<<i<<' '<<j<<' '<<k<<' '<<bx[i][j][k]<<' '<<by[i][j][k]<<' '<<bz[i][j][k]<<endl;

                            uxtemp+=ux[i][j][k]/(nx*ny*nz);
                            uytemp+=uy[i][j][k]/(nx*ny*nz);
                            uztemp+=uz[i][j][k]/(nx*ny*nz);

                            bxtemp+=bx[i][j][k]/(nx*ny*nz);
                            bytemp+=by[i][j][k]/(nx*ny*nz);
                            bztemp+=bz[i][j][k]/(nx*ny*nz);


                            if ((i==nx)&&(j==ny)&&(k==nz)){
                            //cout.precision(6);
                            //cout<<"time, i, j, k, bxtemp, bytemp, bztemp = "<<it<<' '<<i<<' '<<j<<' '<<k<<' '<<uxtemp<<' '<<uytemp<<' ' <<uztemp<<' '<<bxtemp<<' '<<bytemp<<' '<<bztemp<<endl;
                            printf("%4d %5.2lf %4d %4d %4d %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf\n", it, time, i, j, k, rho[i][j][k], p[i][j][k],
                                   uxtemp, uytemp, uztemp, ax[i][j][k], ay[i][j][k], az[i][j][k], bxtemp, bytemp, bztemp);
                            }


                            /*
                            if ((i==nx)&&(j==ny)&&(k==nz)){
                            fp=fopen("data.txt", "a+");
                            fprintf(fp, "\n%5.3lf %4d %4d %4d %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf", it, i, j, k, uxtemp, uytemp, uztemp, bxtemp, bytemp, bztemp);
                            fclose(fp);
                            }

                            */


                            rho[i][j][k]+=dt*(-0.5*rho[i][j][k]*((ux[i+1][j][k]-ux[i-1][j][k])/dx+(uy[i][j+1][k]-uy[i][j-1][k])/dy+(uz[i][j][k+1]-uz[i][j][k-1])/dz)
                            -0.5*(ux[i][j][k]*(rho[i+1][j][k]-rho[i-1][j][k])/dx + uy[i][j][k]*(rho[i][j+1][k]-rho[i][j-1][k])/dy+uz[i][j][k]*(rho[i][j][k+1]-rho[i][j][k-1])/dz));
                            //cout<<"t, i, j, k, rho = "<<t<<' '<<i<<' '<<j<<' '<<k<<' '<<rho[i][j][k]<<endl;

                            p[i][j][k]+=dt*(-0.5*gamma*p[i][j][k]*((ux[i+1][j][k]-ux[i-1][j][k])/dx+(uy[i][j+1][k]-uy[i][j-1][k])/dy+(uz[i][j][k+1]-uz[i][j][k-1])/dz)
                            -0.5*(ux[i][j][k]*(p[i+1][j][k]-p[i-1][j][k])/dx + uy[i][j][k]*(p[i][j+1][k]-p[i][j-1][k])/dy + uz[i][j][k]*(p[i][j][k+1]-p[i][j][k-1])/dz));
                            //cout<<"time, i, j, k, p = "<<it<<' '<<i<<' '<<j<<' '<<k<<' '<<p[i][j][k]<<endl;

                            ux[i][j][k]+=dt*(-0.5*ux[i][j][k]*(ux[i+1][j][k]-ux[i-1][j][k])/dx-0.5*uy[i][j][k]*(ux[i][j+1][k]-ux[i][j-1][k])/dy-0.5*uz[i][j][k]*(ux[i][j][k+1]-ux[i][j][k-1])/dz
                            -0.5*(p[i+1][j][k]-p[i-1][j][k])/(dx*rho[i][j][k])
                            -0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(by[i+1][j][k]-by[i-1][j][k])/dy+bz[i][j][k]*(bz[i+1][j][k]-bz[i-1][j][k])/dz)/(mu*rho[i][j][k])
                            +0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(bx[i][j+1][k]-bx[i][j-1][k])/dy+bz[i][j][k]*(bx[i][j][k+1]-bx[i][j][k-1])/dz)/(mu*rho[i][j][k])
                            +nu*((ux[i+1][j][k]+ux[i-1][j][k]-2.0*ux[i][j][k])/(dx*dx)+(ux[i][j+1][k]+ux[i][j-1][k]-2.0*ux[i][j][k])/(dy*dy)+(ux[i][j][k+1]+ux[i][j][k-1]-2.0*ux[i][j][k])/(dz*dz))/(rho[i][j][k])+fvx);
                            //cout<<"t, i, j, k, ux = "<<t<<' '<<i<<' '<<j<<' '<<k<<' '<<ux[i][j][k]<<endl;

                            uy[i][j][k]+=dt*(-0.5*ux[i][j][k]*(uy[i+1][j][k]-uy[i-1][j][k])/dx-0.5*uy[i][j][k]*(uy[i][j+1][k]-uy[i][j-1][k])/dy-0.5*uz[i][j][k]*(uy[i][j][k+1]-uy[i][j][k-1])/dz
                            -0.5*(p[i][j+1][k]-p[i][j-1][k])/(dy*rho[i][j][k])
                            -0.5*(bx[i][j][k]*(bx[i][j+1][k]-bx[i][j-1][k])/dy+by[i][j][k]*(by[i][j+1][k]-by[i][j-1][k])/dy+bz[i][j][k]*(bz[i][j+1][k]-bz[i][j-1][k])/dz)/(mu*rho[i][j][k])
                            +0.5*(bx[i][j][k]*(by[i+1][j][k]-by[i-1][j][k])/dx+by[i][j][k]*(by[i][j+1][k]-by[i][j-1][k])/dy+bz[i][j][k]*(by[i][j][k+1]-by[i][j][k-1])/dz)/(mu*rho[i][j][k])
                            +nu*((uy[i+1][j][k]+uy[i-1][j][k]-2.0*uy[i][j][k])/(dx*dx)+(uy[i][j+1][k]+uy[i][j-1][k]-2.0*uy[i][j][k])/(dy*dy)+(uy[i][j][k+1]+uy[i][j][k-1]-2.0*uy[i][j][k])/(dz*dz))/(rho[i][j][k])+fvy);
                            //cout<<"time, i, j, k, uy = "<<it<<' '<<i<<' '<<j<<' '<<k<<' '<<uy[i][j][k]<<endl;

                            uz[i][j][k]+=dt*(-0.5*ux[i][j][k]*(uz[i+1][j][k]-uz[i-1][j][k])/dx-0.5*uy[i][j][k]*(uz[i][j+1][k]-uz[i][j-1][k])/dy-0.5*uz[i][j][k]*(uz[i][j][k+1]-uz[i][j][k-1])/dz
                            -0.5*(p[i][j][k+1]-p[i][j][k-1])/(dz*rho[i][j][k])
                            -0.5*(bx[i][j][k]*(bx[i][j][k+1]-bx[i][j][k-1])/dz+by[i][j][k]*(by[i][j][k+1]-by[i][j][k-1])/dz+bz[i][j][k]*(bz[i][j][k+1]-bz[i][j][k-1])/dz)/(mu*rho[i][j][k])
                            +0.5*(bx[i][j][k]*(bz[i+1][j][k]-bz[i-1][j][k])/dx+by[i][j][k]*(bz[i][j+1][k]-bz[i][j-1][k])/dy+bz[i][j][k]*(bz[i][j][k+1]-bz[i][j][k-1])/dz)/(mu*rho[i][j][k])
                            +nu*((uz[i+1][j][k]+uz[i-1][j][k]-2.0*uz[i][j][k])/(dx*dx)+(uz[i][j+1][k]+uz[i][j-1][k]-2.0*uz[i][j][k])/(dy*dy)+(uz[i][j][k+1]+uz[i][j][k-1]-2.0*uz[i][j][k])/(dz*dz))/(rho[i][j][k])+fvz);
                            //cout<<"t, i, j, k, uz = "<<t<<' '<<i<<' '<<j<<' '<<k<<' '<<uz[i][j][k]<<endl;

                            ax[i][j][k]+=dt*(0.5*uy[i][j][k]*((ay[i+1][j][k]-ay[i-1][j][k])/dx-(ax[i][j+1][k]-ax[i][j-1][k])/dy)
                            -0.5*uz[i][j][k]*((ax[i][j][k+1]-ax[i][j][k-1])/dz-(az[i+1][j][k]-az[i-1][j][k])/dx)
                            +eta*((ax[i+1][j][k]+ax[i-1][j][k]-2.0*ax[i][j][k])/(dx*dx)+(ax[i][j+1][k]+ax[i][j-1][k]-2.0*ax[i][j][k])/(dy*dy)+(ax[i][j][k+1]+ax[i][j][k-1]-2.0*ax[i][j][k])/(dz*dz)));
                            //cout<<"t, i, j, k, ax = "<<t<<' '<<i<<' '<<j<<' '<<k<<' '<<ax[i][j][k]<<endl;

                            ay[i][j][k]+=dt*(0.5*uz[i][j][k]*((az[i][j+1][k]-az[i][j-1][k])/dy-(ay[i][j][k+1]-ay[i][j][k-1])/dz)
                            -0.5*ux[i][j][k]*((ay[i+1][j][k]-ay[i-1][j][k])/dx-(ax[i][j+1][k]-ax[i][j-1][k])/dy)
                            +eta*((ay[i+1][j][k]+ay[i-1][j][k]-2.0*ay[i][j][k])/(dx*dx)+(ay[i][j+1][k]+ay[i][j-1][k]-2.0*ay[i][j][k])/(dy*dy)+(ay[i][j][k+1]+ay[i][j][k-1]-2.0*ay[i][j][k])/(dz*dz)));
                            //cout<<"t, i, j, k, ay = "<<t<<' '<<i<<' '<<j<<' '<<k<<' '<<ay[i][j][k]<<endl;

                            az[i][j][k]+=dt*(0.5*ux[i][j][k]*((ax[i][j][k+1]-ax[i][j][k-1])/dz-(az[i+1][j][k]-az[i-1][j][k])/dx)
                            -0.5*uy[i][j][k]*((az[i][j+1][k]-az[i][j-1][k])/dy-(ay[i][j][k+1]-ay[i][j][k-1])/dz)
                            +eta*((az[i+1][j][k]+az[i-1][j][k]-2.0*az[i][j][k])/(dx*dx)+(az[i][j+1][k]+az[i][j-1][k]-2.0*az[i][j][k])/(dy*dy)+(az[i][j][k+1]+az[i][j][k-1]-2.0*az[i][j][k])/(dz*dz)));
                            //cout<<"t, i, j, k, az = "<<it<<' '<<i<<' '<<j<<' '<<k<<' '<<az[i][j][k]<<endl;

                     }
              }
       }

       time+=dt;

       //cout<<'test bx = '<<bx[2][1][1]<<endl;

}


 return 0;

}
