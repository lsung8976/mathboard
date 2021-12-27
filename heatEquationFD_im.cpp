//유한차분법 적용코드//함축적 방법

????????
#include <iostream>
#include <math.h>
#include <fstream>
#define _USE_MATH_DEFINES
using namespace std;

int main()
{
    // 파일 쓰기 준비
    std::ofstream out("test3.txt", std::ios::app);
    
    //열방정식 구하기

    //권리행사가격 : E = 230
    double alpha = 2;
    double xmult;
    //시간 : t=0(현재시점) ~ T=0.125(만기시점)
    double T = 0.1;

    //격자간격 Nx = 30;
    int Nx = 12;
    //시간간격 Nt = round(T/k);
    int Nt = 5;
    //증분
    int h;

    //격자 배열:
    double x_arr[Nx+2];
    
    double dd[Nx+2], c[Nx+2], a[Nx+2], d[Nx+2], b[Nx+2];
    double u_arr[Nt+2][Nx+2];

    
    //u배열 초기화 (경계조건 겸용 초기화)
    for(int i = 1;i <= Nt + 1; i++)
    {
        for(int j = 1;j <= Nx + 1; j++)
        {
            u_arr[i][j] = 0.0;
        }
    }
    
    //간격 크기
    double size_x = (1.0) / Nx;
    double size_t = T / Nt;
    //x 초기화
    for(int i = 1; i <= Nx - 2; i++)
    {
        x_arr[i] = i * size_x;
        dd[i] = 1 + 2 * alpha;
        c[i] = - alpha;
        a[i] = - alpha;
        b[i] = 0;
    }
    
    //증분 설정 // x(2) - x(1)
    h = x_arr[2] - x_arr[1];
    //k
    double k = alpha * h * h;
    
    //초기값 설정 (초기조건) 
    for(int i = 1; i <= Nx + 1; i++)   
        u_arr[1][i] = sin(M_PI * x_arr[i]);
    

    
    for(int n = 1; n <= Nt; n++)
    {
        for(int i = 1; i <= Nx - 2; i++)
        {
            d[i] = dd[i];
            if(i != 1)
            {
                xmult = a[i-1]/d[i-1];
                d[i] = d[i] - xmult * c[i - 1];
                b[i] = b[i] - xmult * b[i - 1];
            }
        }

        u_arr[n+1][Nx-1] = b[Nx-2] / d[Nx-2];

        for(int i = 1 ; i <= Nt + 1; i++)   //시간
        {
            for(int j = 1; j <= Nx + 1; j++)
            {
                cout << u_arr[i][j] << " ";
            }
            cout << endl;
        }
        // for(int i = Nx - 3; i >= 0 ; i--)
        // {
        //     u_arr[n+1][i+1] = (b[i] - c[i]*u_arr[n+1][i+2])/d[i];
        //     cout << u_arr[n+1][i+2] << " ";
        // }
        // cout << endl;
    }
    
    
    // for(int i = 1 ; i <= Nt + 1; i++)   //시간
    // {
    //     for(int j = 1; j <= Nx + 1; j++)
    //     {
    //         cout << u_arr[i][j] << " ";
    //     }
       
    //     cout << endl;
    // }

    // for(int i = 1 ; i <= Nt + 1; i++)   //시간
    // {
    //     for(int j = 1; j <= Nx + 1; j++)
    //     {
    //         if (out.is_open()) {
    //             out << i * size_t << " " << x_arr[j] << " " << u_arr[i][j] << endl;
    //         }
    //     }
    // }
    
    return 0;
}