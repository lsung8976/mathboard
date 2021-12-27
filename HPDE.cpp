// Hyperbolic Partial differential Equation 쌍곡형 미분방정식
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct flow{
    double old;
    double current;
    double newData;
};

double f(double x){
    if((x >= 0.0) && (x <= 1.0))
        return 0.005*x;
    else
        return -0.01*x + 0.015;
}

int main(){
    int i, count;
    int print_step, mesh;
    struct flow displace[30];
    double x;
    double dx, dt, sqr_c, r;
    double time, end_time;

    printf("\nWrite x range from 0 to L & mesh Size : ");
    scanf("%lf%d", &x, &mesh);
    printf("Input the differential equation coefficient : ");
    scanf("%lf", &sqr_c);
    printf("Input elapsed time to end caculation : ");
    scanf("%lf", &end_time);
    printf("Input iteration numbers to print out the results : ");
    scanf("%d", &print_step);

    dx = x / mesh;
    dt = dx / sqrt(sqr_c);
    printf("\n The time step(dt) is %lf \n Input the new dt : ", dt);
    scanf("%lf", &dt);
    r = sqr_c * dt * dt / (dx * dx);
    time = dt;
    count = 1;

    printf("\n time");
    for (i=0;i<mesh;i+=2)
        printf("    %2d", i + 1);
    if(!(i%2)) printf("     %2d", mesh + 1);
    printf("\n\n");

    displace[0].old = displace[0].current = displace[0].newData = 0.0;
    displace[mesh] = displace[0];

    for(i = 0; i <= mesh; i++)
        displace[i].current = f(i*dx);  //u(i,j) = f(i * dx)
    for(i = 1; i < mesh; i++)   //u(i,j+1) = (1-r)u(i,j) + 1/2 * (u(i+1) + u(i-1))
        displace[i].newData = (1.0 - r) * displace[i].current + 0.5 * r * (displace[i+1].current+displace[i-1].current);
    do{
        for(i = 1; i < mesh; i++){
            displace[i].old = displace[i].current;
            displace[i].current = displace[i].newData;  //새 값이 들어왔으니 현재 값 -> 올드값으로 바뀌고 새 값 -> 현재 값으로 바뀜
        }
        for(i = 1; i < mesh; i++){
            displace[i].newData = 2.0*(1.0 - r) * displace[i].current + r * (displace[i - 1].current + displace[i+1].current) - displace[i].old;
        }
        time += dt;
        if((++count%print_step) == 0){
            printf("%7.5lf", time);
            for(i = 0; i < mesh; i++)
                printf("%8.4lf", displace[i].newData);
            if(!(i%2))
                printf("%8.4lf", displace[mesh].newData);
        }
    }while((time - print_step) != 0);
    
    if((count % print_step) != 0){
        printf("%7.5lf", time);
        for(i = 0 ;i < mesh; i++)
            printf("%8.4lf", displace[i].newData);
        if(!(i%2))
                printf("%8.4lf", displace[mesh].newData);
    }

    return 0;
}
