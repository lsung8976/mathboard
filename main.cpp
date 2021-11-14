#include<iostream>

#define nx 10
#define ny 10
#define nz 10
#define gamma 1.667
#define dt 1.0
#define dx 1.0
#define dy 1.0
#define dz 1.0

using namespace std;

int main()
{

int i, j, k, n, nstep=30;

double rho[nx+2][ny+2][nz+2]={0.0}, p[nx+2][ny+2][nz+2]={0.0};
double ux[nx+2][ny+2][nz+2]={0.0}, uy[nx+2][ny+2][nz+2]={0.0}, uz[nx+2][ny+2][nz+2]={0.0};
double ax[nx+2][ny+2][nz+2]={0.0}, ay[nx+2][ny+2][nz+2]={0.0}, az[nx+2][ny+2][nz+2]={0.0};
double bx[nx+2][ny+2][nz+2]={0.0}, by[nx+2][ny+2][nz+2]={0.0}, bz[nx+2][ny+2][nz+2]={0.0};

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

double mu=1.0, nu=1.0, eta=1.0;



k_rho[2][2][2]=10.2;
k_ux[1][1][1]=10.3;
k_az[1][1][1]=10.3;

//cout<<"No  x     y[0]      y[1]"<<endl;

for(n=1; n<=nstep; ++n){

	/*

	for(i=1; i<=nx; ++i){
	for(j=1; j<=ny; ++j){
	for(k=1; k<=nz; ++k){
	ax[1][j][k]=0.345;
	ay[i][2][k]=0.2;
	az[i][j][1]=0.2;
	}
	}
	}

	*/

	for(i=1; i<=nx; ++i){
		for(j=1; j<=ny; ++j){
			for(k=1; k<=nz; ++k){
			bx[i][j][k]=0.5*(az[i][j+1][k]-az[i][j-1][k])/dy - 0.5*(ay[i][j][k+1]-ay[i][j][k-1])/dz;
			by[i][j][k]=0.5*(ax[i][j][k+1]-ax[i][j][k-1])/dz - 0.5*(az[i+1][j][k]-az[i-1][j][k])/dx;
			bz[i][j][k]=0.5*(ay[i+1][j][k]-ay[i-1][j][k])/dx - 0.5*(ax[i][j+1][k]-ax[i][j-1][k])/dy;
			}
		}
	}


	for(i=1; i<=nx; ++i){
		for(j=1; j<=ny; ++j){
			for(k=1; k<=nz; ++k){
			k_rho[i][j][k]=-0.5*rho[i][j][k]*((ux[i+1][j][k]-ux[i-1][j][k])/dx+(uy[i][j+1][k]-uy[i][j-1][k])/dy+(uz[i][j][k+1]-uz[i][j][k-1])/dz)
			-0.5*(ux[i][j][k]*(rho[i+1][j][k]-rho[i-1][j][k])/dx + uy[i][j][k]*(rho[i][j+1][k]-rho[i][j-1][k])/dy+uz[i][j][k]*(rho[i][j][k+1]-rho[i][j][k-1])/dz);
			//cout<<i<<' '<<j<<' '<<k<<' '<<rho[i][j][k]<<endl;
			}
		}
	}


	for(i=1; i<=nx; ++i){
		for(j=1; j<=ny; ++j){
			for(k=1; k<=nz; ++k){
			k_p[i][j][k]=-0.5*gamma*p[i][j][k]*((ux[i+1][j][k]-ux[i-1][j][k])/dx+(uy[i][j+1][k]-uy[i][j-1][k])/dy+(uz[i][j][k+1]-uz[i][j][k-1])/dz)
			-0.5*(ux[i][j][k]*(p[i+1][j][k]-p[i-1][j][k])/dx + uy[i][j][k]*(p[i][j+1][k]-p[i][j-1][k])/dy + uz[i][j][k]*(p[i][j][k+1]-p[i][j][k-1])/dz);
			}
		}
	}


	for(i=1; i<=nx; ++i){
		for(j=1; j<=ny; ++j){
			for(k=1; k<=nz; ++k){
			k_ux[i][j][k]=-0.5*ux[i][j][k]*(ux[i+1][j][k]-ux[i-1][j][k])/dx-0.5*uy[i][j][k]*(ux[i][j+1][k]-ux[i][j-1][k])/dy-0.5*uz[i][j][k]*(ux[i][j][k+1]-ux[i][j][k-1])/dz
			-0.5*(p[i+1][j][k]-p[i-1][j][k])/(dx*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(by[i+1][j][k]-by[i-1][j][k])/dy+bz[i][j][k]*(bz[i+1][j][k]-bz[i-1][j][k])/dz)/(mu*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(bx[i][j+1][k]-bx[i][j-1][k])/dy+bz[i][j][k]*(bx[i][j][k+1]-bx[i][j][k-1])/dz)/(mu*rho[i][j][k])
			+nu*((ux[i+1][j][k]+ux[i-1][j][k]-2.0*ux[i][j][k])/(dx*dx)+(ux[i][j+1][k]+ux[i][j-1][k]-2.0*ux[i][j][k])/(dy*dy)+(ux[i][j][k+1]+ux[i][j][k-1]-2.0*ux[i][j][k])/(dz*dz))/(rho[i][j][k]);

			/*
			k_ux[i][j][k]=-0.5*ux[i][j][k]*(ux[i+1][j][k]-ux[i-1][j][k])/dx-0.5*uy[i][j][k]*(ux[i][j+1][k]-ux[i][j-1][k])/dy-0.5*uz[i][j][k]*(ux[i][j][k+1]-ux[i][j][k-1])/dz
			-0.5*(p[i+1][j][k]-p[i-1][j][k])/(dx*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(by[i+1][j][k]-by[i-1][j][k])/dy+bz[i][j][k]*(bz[i+1][j][k]-bz[i-1][j][k])/dz)/(mu*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i+1][j][k]-bx[i-1][j][k])/dx+by[i][j][k]*(bx[i][j+1][k]-bx[i][j-1][k])/dy+bz[i][j][k]*(bx[i][j][k+1]-bx[i][j][k-1])/dz)/(mu*rho[i][j][k])
			+nu*((ux[i+1][j][k]+ux[i-1][j][k]-2.0*ux[i][j][k])/(dx*dx)+(ux[i][j+1][k]+ux[i][j-1][k]-2.0*ux[i][j][k])/(dy*dy)+(ux[i][j][k+1]+ux[i][j][k-1]-2.0*ux[i][j][k])/(dz*dz))/(rho[i][j][k]);
			*/

			k_uy[i][j][k]=-0.5*ux[i][j][k]*(uy[i+1][j][k]-uy[i-1][j][k])/dx-0.5*uy[i][j][k]*(uy[i][j+1][k]-uy[i][j-1][k])/dy-0.5*uz[i][j][k]*(uy[i][j][k+1]-uy[i][j][k-1])/dz
			-0.5*(p[i][j+1][k]-p[i][j-1][k])/(dy*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i][j+1][k]-bx[i][j-1][k])/dy+by[i][j][k]*(by[i][j+1][k]-by[i][j-1][k])/dy+bz[i][j][k]*(bz[i][j+1][k]-bz[i][j-1][k])/dz)/(mu*rho[i][j][k])
			+0.5*(bx[i][j][k]*(by[i+1][j][k]-by[i-1][j][k])/dx+by[i][j][k]*(by[i][j+1][k]-by[i][j-1][k])/dy+bz[i][j][k]*(by[i][j][k+1]-by[i][j][k-1])/dz)/(mu*rho[i][j][k])
			+nu*((uy[i+1][j][k]+uy[i-1][j][k]-2.0*uy[i][j][k])/(dx*dx)+(uy[i][j+1][k]+uy[i][j-1][k]-2.0*uy[i][j][k])/(dy*dy)+(uy[i][j][k+1]+uy[i][j][k-1]-2.0*uy[i][j][k])/(dz*dz))/(rho[i][j][k]);

			k_uz[i][j][k]=-0.5*ux[i][j][k]*(uz[i+1][j][k]-uz[i-1][j][k])/dx-0.5*uy[i][j][k]*(uz[i][j+1][k]-uz[i][j-1][k])/dy-0.5*uz[i][j][k]*(uz[i][j][k+1]-uz[i][j][k-1])/dz
			-0.5*(p[i][j][k+1]-p[i][j][k-1])/(dz*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bx[i][j][k+1]-bx[i][j][k-1])/dz+by[i][j][k]*(by[i][j][k+1]-by[i][j][k-1])/dz+bz[i][j][k]*(bz[i][j][k+1]-bz[i][j][k-1])/dz)/(mu*rho[i][j][k])
			+0.5*(bx[i][j][k]*(bz[i+1][j][k]-bz[i-1][j][k])/dx+by[i][j][k]*(bz[i][j+1][k]-bz[i][j-1][k])/dy+bz[i][j][k]*(bz[i][j][k+1]-bz[i][j][k-1])/dz)/(mu*rho[i][j][k])
			+nu*((uz[i+1][j][k]+uz[i-1][j][k]-2.0*uz[i][j][k])/(dx*dx)+(uz[i][j+1][k]+uz[i][j-1][k]-2.0*uz[i][j][k])/(dy*dy)+(uz[i][j][k+1]+uz[i][j][k-1]-2.0*uz[i][j][k])/(dz*dz))/(rho[i][j][k]);
			}
		}
	}


	for(i=1; i<=nx; ++i){
		for(j=1; j<=ny; ++j){
			for(k=1; k<=nz; ++k){

			k_ax[i][j][k]=0.5*uy[i][j][k]*((ay[i+1][j][k]-ay[i-1][j][k])/dx-(ax[i][j+1][k]-ax[i][j-1][k])/dy)
			-0.5*uz[i][j][k]*((ax[i][j][k+1]-ax[i][j][k-1])/dz-(az[i+1][j][k]-az[i-1][j][k])/dx)
			+eta*((ax[i+1][j][k]+ax[i-1][j][k]-2.0*ax[i][j][k])/(dx*dx)+(ax[i][j+1][k]+ax[i][j-1][k]-2.0*ax[i][j][k])/(dy*dy)+(ax[i][j][k+1]+ax[i][j][k-1]-2.0*ax[i][j][k])/(dz*dz));

			k_ay[i][j][k]=0.5*uz[i][j][k]*((az[i][j+1][k]-az[i][j-1][k])/dy-(ay[i][j][k+1]-ay[i][j][k-1])/dz)
			-0.5*ux[i][j][k]*((ay[i+1][j][k]-ay[i-1][j][k])/dx-(ax[i][j+1][k]-ax[i][j-1][k])/dy)
			+eta*((ay[i+1][j][k]+ay[i-1][j][k]-2.0*ay[i][j][k])/(dx*dx)+(ay[i][j+1][k]+ay[i][j-1][k]-2.0*ay[i][j][k])/(dy*dy)+(ay[i][j][k+1]+ay[i][j][k-1]-2.0*ay[i][j][k])/(dz*dz));

			k_az[i][j][k]=0.5*ux[i][j][k]*((ax[i][j][k+1]-ax[i][j][k-1])/dz-(az[i+1][j][k]-az[i-1][j][k])/dx)
			-0.5*uy[i][j][k]*((az[i][j+1][k]-az[i][j-1][k])/dy-(ay[i][j][k+1]-ay[i][j][k-1])/dz)
			+eta*((az[i+1][j][k]+az[i-1][j][k]-2.0*az[i][j][k])/(dx*dx)+(az[i][j+1][k]+az[i][j-1][k]-2.0*az[i][j][k])/(dy*dy)+(az[i][j][k+1]+az[i][j][k-1]-2.0*az[i][j][k])/(dz*dz));

			//cout<<"k_az = "<<' '<<k_az[i][j][k]<<' '<<endl;
			}
		}
	}








	//cout<<'test bx = '<<bx[2][1][1]<<endl;

}


 return 0;

}

