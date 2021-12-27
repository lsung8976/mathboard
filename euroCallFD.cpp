//유한차분법 적용코드//명시적 방법

#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;

int main()
{

    // 파일 쓰기 준비
    std::ofstream out("test.txt", std::ios::app);
    
    //유러피안 콜옵션 가격 구하기

    //(고려X)현재의 주가지수 : S = 240
    //double S = 240;
    //권리행사가격 : E = 230
    double E = 230;
    //옵션의 기간 : t=0(현재시점) ~ T=1(만기시점)
    //double T = 2.0 / 12.0;
    int T = 1;
    //주가변동 : sigma = 0.5 //5%
    double sigma = 0.5;
    //비위험이자율 : r = 0.03 //3%
    double r = 0.03;
    //격자길이 0~800
    int L = 800;
    //격자간격 Nx = 50;
    int Nx = 50;
    //시간간격 Nt = 1000
    int Nt = 1000;
    //증분
    double h;

    //k
    double k = double(T) / Nt;
    double t_arr[Nt+2];

    //격자 배열:
    double x_arr[Nx+2];
    double u_arr[Nt+2][Nx+2];

    
    //u배열 초기화
    for(int i = 1;i <= Nt + 1; i++)
    {
        for(int j = 1;j <= Nx + 1; j++)
        {
            u_arr[i][j] = 0;
        }
    }
    
    
    //간격 크기:
    int size_x = L / Nx;
    for(int i = 1; i <= Nx; i++)
    {
        x_arr[i] = i * size_x;
    }
    
    //증분 설정 // x(2) - x(1)
    h = x_arr[2] - x_arr[1];

    for(int i = 1; i <= Nx; i++)    //초기값 설정
    {
        if (x_arr[i] <= E) 
        {   
            u_arr[1][i] = 0; 
        }
        else {
            u_arr[1][i] = x_arr[i] - E;
        }
    }
    
    for(int i = 2; i <= Nt + 1; i++)    //초기값 설정
        u_arr[i][Nx] = L - E * exp(-r * k * (i - 1));
    

    
    for(int n = 1; n <= Nt; n++)
    {
        for(int i = 2 ; i <= Nx - 1;i++)
        {
            u_arr[n+1][i] = u_arr[n][i] + k * ((1/2) * (sigma*sigma) * (double(i-1)*h)*(double(i-1)*h) 
                * ((u_arr[n][i+1] - 2 *u_arr[n][i] + u_arr[n][i-1]) / (h*h)) 
                + r * double(i-1) * h * ((u_arr[n][i+1] - u_arr[n][i-1])/(2*h)) - r * u_arr[n][i] );
        }
    }
    for(int i = 1 ; i <= Nt + 1; i++)
        t_arr[i] = double(i) * k;

    
    for(int i = 1 ; i <= Nt + 1; i++)
    {
        for(int j = 1; j <= Nx + 1; j++)
        {
            if (out.is_open()) {
                out << t_arr[i] << " " << x_arr[j] << " " << u_arr[i][j] << endl;
            }
        }
    }
    
    return 0;
}