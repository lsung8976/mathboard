#include <stdio.h>
#include <math.h>

//누적정규분포 (C++ 에는 누적표준 정규분포가 없음 -> 근사값인 0.5(1-erfc)를 활용함 (erf : error function))
double normalCDF(double value)
{
        return 0.5 * erfc(-value * M_SQRT1_2);
}

// https://www.johndcook.com/blog/cpp_phi/ 이곳에서 참고한 누적정규분포 근사함수
double phi(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

int main()
{
    double d1, temp_d1, d2, CallPrice1, CallPrice2;
    //유러피안 콜옵션 가격 구하기
    //현재의 주가지수 : S = 240
    double S = 240;
    //권리행사가격 : E = 250
    double E = 250;
    //옵션의 기간 : T = 2/12 // 2개월
    double T = 2.0 / 12.0;
    //주가변동 : sigma = 0.38 //38%
    double sigma = 0.38;
    //비위험이자율 : r = 0.06 //6%
    double r = 0.06;

    
    temp_d1 = (log10(S/E) + (r + 0.5 * (sigma * sigma) * T));  // /(sigma*sqrt(T));
    d1 = temp_d1 / (sigma*sqrt(T));
    d2 = d1 - (sigma*sqrt(T));

    CallPrice1 = S * normalCDF(d1) - E * exp(-r * T)*normalCDF(d2);
    CallPrice2 = S * phi(d1) - E * exp(-r * T)*phi(d2);

    printf("CallPrice1 : %lf\n", CallPrice1);
    printf("CallPrice2 : %lf\n", CallPrice2);

    //printf("%lf\n",log10(S/E));
    //printf("%lf\n", (sigma*sqrt(T)));

    printf("%lf\n", normalCDF(2));
    printf("%lf\n", phi(2));

    //printf("d1 d2 : %lf %lf", d1, d2);
    return 0;
}