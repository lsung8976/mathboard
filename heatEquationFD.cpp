//유한차분법 적용코드//명시적 방법

#include <iostream>
#include <math.h>
#include <fstream>
#define _USE_MATH_DEFINES
using namespace std;

int main()
{
    // 파일 쓰기 준비
    std::ofstream out("test2.txt", std::ios::app);
    
    //열방정식 구하기

    //권리행사가격 : E = 230
    double alpha = 0.45;
    //시간 : t=0(현재시점) ~ T=0.125(만기시점)
    double T = 0.125;

    

    //k
    double k = pow(alpha * 2, 2);

    //격자간격 Nx = 30;
    int Nx = 30;
    //시간간격 Nt = round(T/k);
    int Nt = 500;
    //증분
    int h;

    //격자 배열:
    double x_arr[Nx+2];
    double u_arr[Nt+2][Nx+2], ex_u_arr[Nt+2][Nx+2];;

    
    //u배열 초기화 (경계조건 겸용 초기화)
    for(int i = 1;i <= Nt + 1; i++)
    {
        for(int j = 1;j <= Nx + 1; j++)
        {
            u_arr[i][j] = 0.0;
        }
    }
    
    //간격 크기:
    double size_x = (1.0) / Nx;
    double size_t = T / Nt;
    for(int i = 1; i <= Nx + 1; i++)
    {
        x_arr[i] = i * size_x;
        //cout << x_arr[i] << " ";
    }
    //cout << endl;
    
    //증분 설정 // x(2) - x(1)
    h = x_arr[2] - x_arr[1];
    
    for(int i = 1; i <= Nx + 1; i++)    //초기값 설정 (초기조건) 
        u_arr[1][i] = sin(M_PI * x_arr[i]);
    

    
    for(int n = 1; n <= Nt; n++)
    {
        for(int i = 2 ; i <= Nx - 1;i++)
        {
            u_arr[n+1][i] = u_arr[n][i] + alpha * ( (u_arr[n][i-1] - 2*u_arr[n][i] + u_arr[n][i+1]) );
            ex_u_arr[n+1][i] = sin(M_PI * x_arr[i]) * exp( -(M_PI*M_PI) * (k * n) );
        }
    }


    
    // for(int i = 1 ; i <= Nt + 1; i++)   //시간
    // {
    //     for(int j = 1; j <= Nx + 1; j++)
    //     {
    //         cout << u_arr[i][j] << " ";
    //     }
       
    //     cout <<endl;
    // }

    for(int i = 1 ; i <= Nt + 1; i++)   //시간
    {
        for(int j = 1; j <= Nx + 1; j++)
        {
            if (out.is_open()) {
                out << i * size_t << " " << x_arr[j] << " " << u_arr[i][j] << endl;
            }
        }
    }
    
    return 0;
}