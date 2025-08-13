#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>

#define  GAMMA                1.67
#define  max_no_variable      8
#define  L                    0.8
#define  mu                   1
#define  ghost                6
#define  cells                450
#define  faces                cells
#define  points               cells+faces+ghost+7
#define  size                 points-1


#define  cr1         0.3333
#define  cr2         1.16666
#define  cr3         1.8333
#define  cr4          .83333
#define   cr5         .166666
 
#define   gamma1     0.1
#define   gamma2     0.6
#define   gamma3     0.3
 
#define    elipson    1e-06

 
using namespace std;


int i,j,k;
 
double  SI1,SI2,SI3,SI4,omega1,omega2,omega3,omega_1,omega_2,omega_3,w1,w2,w3,w_1,w_2,w_3,S1,S2,S3,S_1,S_2,S_3;

double pressure_L[points][points],pressure_R[points][points],ethalpy_R[points][points],density_L[points][points],density_R[points][points],ethalpy_L[points][points],velocity_yR[points][points],velocity_yL[points][points];                               

double velocity_xR[points][points],velocity_xL[points][points],magnetic_L[points][points],magnetic_R[points][points],bx_R[points][points],bx_L[points][points];

double p_L[points][points],p_R[points][points],e_R[points][points],d_L[points][points],d_R[points][points],e_L[points][points],v_yR[points][points],v_yL[points][points];                               

double v_xR[points][points],v_xL[points][points],BY_L[points][points],BY_R[points][points],BX_R[points][points],BX_L[points][points], bz_R[points][points],bz_L[points][points], BZ_R[points][points], BZ_L[points][points] ;

double velocity_zL[points][points], velocity_zR[points][points] , v_zL[points][points], v_zR[points][points];


void WENO1(double u[][points]);
void WENO2(double u1[][points]);
void WENO3(double u2[][points]);
void WENO4(double u3[][points]);
void WENO5(double u4[][points]);
void WENO6(double u5[][points]); 
void WENO7(double u6[][points]); 
void WENO8(double u7[][points]);  
void WENO9(double u8[][points]);  
 


 void WENO1(double u[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for x-direction
 
 SI1 = (1.0833)*(u[i][j-4]-2*u[i][j-2]+u[i][j])*(u[i][j-4]-2*u[i][j-2]+u[i][j]) + 0.25*(u[i][j-4]-4*u[i][j-2]+3*u[i][j])*(u[i][j-4]-4*u[i][j-2]+3*u[i][j]);
 
 SI2 = (1.08333)*(u[i][j-2]-2*u[i][j]+u[i][j+2])*(u[i][j-2]-2*u[i][j]+u[i][j+2]) + 0.25*(u[i][j-2]-u[i][j+2])*(u[i][j-2]-u[i][j+2]);
 
 SI3 = (1.08333)*(u[i][j]-2*u[i][j+2]+u[i][j+4])*(u[i][j]-2*u[i][j+2]+u[i][j+4]) + 0.25*(3*u[i][j]-4*u[i][j+2]+u[i][j+4])*(3*u[i][j]-4*u[i][j+2]+u[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u[i][j-4]- cr2*u[i][j-2] + cr3*u[i][j]; 
 	
 S2 = - cr5*u[i][j-2] + cr4*u[i][j] + cr1*u[i][j+2]; 
 
 S3 = - cr5*u[i][j+4] + cr4*u[i][j+2] + cr1*u[i][j]; 
 	
 	
 S_3 = cr1*u[i][j+4]- cr2*u[i][j+2] + cr3*u[i][j];
 
 S_2 = - cr5*u[i][j+2] + cr4*u[i][j] + cr1*u[i][j-2]; 
 
 S_1 = - cr5*u[i][j-4] + cr4*u[i][j-2] + cr1*u[i][j]; 
 
 
  density_L[i][j+1]=(S1*w1)+(S2*w2)+(S3*w3);	
 	
  density_R[i][j-1]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction reconstruction
  
  
  SI1 = (1.08)*(u[i-4][j]-2*u[i-2][j]+u[i][j])*(u[i-4][j]-2*u[i-2][j]+u[i][j]) + 0.25*(u[i-4][j]-4*u[i-2][j]+3*u[i][j])*(u[i-4][j]-4*u[i-2][j]+3*u[i][j]);
 
 SI2 = (1.08)*(u[i-2][j]-2*u[i][j]+u[i+2][j])*(u[i-2][j]-2*u[i][j]+u[i+2][j]) + 0.25*(u[i-2][j]-u[i+2][j])*(u[i-2][j]-u[i+2][j]);
 
 SI3 = (1.08)*(u[i][j]-2*u[i+2][j]+u[i+4][j])*(u[i][j]-2*u[i+2][j]+u[i+4][j]) + 0.25*(3*u[i][j]-4*u[i+2][j]+u[i+4][j])*(3*u[i][j]-4*u[i+2][j]+u[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u[i-4][j]- cr2*u[i-2][j] + cr3*u[i][j]; 
 	
 S2 = - cr5*u[i-2][j] + cr4*u[i][j] + cr1*u[i+2][j]; 
 
 S3 = - cr5*u[i+4][j] + cr4*u[i+2][j] + cr1*u[i][j]; 
 	
 	
 S_3 = cr1*u[i+4][j]- cr2*u[i+2][j] + cr3*u[i][j];
 
 S_2 = - cr5*u[i+2][j] + cr4*u[i][j] + cr1*u[i-2][j]; 
 
 S_1 = - cr5*u[i-4][j] + cr4*u[i-2][j] + cr1*u[i][j]; 
 
  
  
  d_L[i+1][j]=(S1*w1)+(S2*w2)+(S3*w3);	
 	
  d_R[i-1][j]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
  
 	 
	  //if(i==50)
	
 	 }
}

return;
}



 void WENO2(double u1[][points])
 {


for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
 SI1 = (1.0833)*(u1[i][j-4]-2*u1[i][j-2]+u1[i][j])*(u1[i][j-4]-2*u1[i][j-2]+u1[i][j]) + 0.25*(u1[i][j-4]-4*u1[i][j-2]+3*u1[i][j])*(u1[i][j-4]-4*u1[i][j-2]+3*u1[i][j]);
 
 SI2 = (1.08333)*(u1[i][j-2]-2*u1[i][j]+u1[i][j+2])*(u1[i][j-2]-2*u1[i][j]+u1[i][j+2]) + 0.25*(u1[i][j-2]-u1[i][j+2])*(u1[i][j-2]-u1[i][j+2]);
 
 SI3 = (1.08333)*(u1[i][j]-2*u1[i][j+2]+u1[i][j+4])*(u1[i][j]-2*u1[i][j+2]+u1[i][j+4]) + 0.25*(3*u1[i][j]-4*u1[i][j+2]+u1[i][j+4])*(3*u1[i][j]-4*u1[i][j+2]+u1[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u1[i][j-4]- cr2*u1[i][j-2] + cr3*u1[i][j]; 
 	
 S2 = - cr5*u1[i][j-2] + cr4*u1[i][j] + cr1*u1[i][j+2]; 
 
 S3 = - cr5*u1[i][j+4] + cr4*u1[i][j+2] + cr1*u1[i][j]; 
 	
 	
 S_3 = cr1*u1[i][j+4]- cr2*u1[i][j+2] + cr3*u1[i][j];
 
 S_2 = - cr5*u1[i][j+2] + cr4*u1[i][j] + cr1*u1[i][j-2]; 
 
 S_1 = - cr5*u1[i][j-4] + cr4*u1[i][j-2] + cr1*u1[i][j]; 
 
 
  pressure_L[i][j+1]=S1*w1+S2*w2+S3*w3;	
 	
  pressure_R[i][j-1]=S_1*w_1+S_2*w_2+S_3*w_3;	
 	
 	
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u1[i-4][j]-2*u1[i-2][j]+u1[i][j])*(u1[i-4][j]-2*u1[i-2][j]+u1[i][j]) + 0.25*(u1[i-4][j]-4*u1[i-2][j]+3*u1[i][j])*(u1[i-4][j]-4*u1[i-2][j]+3*u1[i][j]);
 
 SI2 = (1.08)*(u1[i-2][j]-2*u1[i][j]+u1[i+2][j])*(u1[i-2][j]-2*u1[i][j]+u1[i+2][j]) + 0.25*(u1[i-2][j]-u1[i+2][j])*(u1[i-2][j]-u1[i+2][j]);
 
 SI3 = (1.08)*(u1[i][j]-2*u1[i+2][j]+u1[i+4][j])*(u1[i][j]-2*u1[i+2][j]+u1[i+4][j]) + 0.25*(3*u1[i][j]-4*u1[i+2][j]+u1[i+4][j])*(3*u1[i][j]-4*u1[i+2][j]+u1[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u1[i-4][j]- cr2*u1[i-2][j] + cr3*u1[i][j]; 
 	
 S2 = - cr5*u1[i-2][j] + cr4*u1[i][j] + cr1*u1[i+2][j]; 
 
 S3 = - cr5*u1[i+4][j] + cr4*u1[i+2][j] + cr1*u1[i][j]; 
 	
 	
 S_3 = cr1*u1[i+4][j]- cr2*u1[i+2][j] + cr3*u1[i][j];
 
 S_2 = - cr5*u1[i+2][j] + cr4*u1[i][j] + cr1*u1[i-2][j]; 
 
 S_1 = - cr5*u1[i-4][j] + cr4*u1[i-2][j] + cr1*u1[i][j]; 
 
  
 	
 p_L[i+1][j]=S1*w1+S2*w2+S3*w3;	
 	
  p_R[i-1][j]=S_1*w_1+S_2*w_2+S_3*w_3;	
 	
 	
 	  //cout<< p_L[i+1][j]<<"\t"<<i<<endl;
 	
	  } 
  
  //return;
  }	 
}

void WENO3(double u2[][points])
{

for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 SI1 = (1.0833)*(u2[i][j-4]-2*u2[i][j-2]+u2[i][j])*(u2[i][j-4]-2*u2[i][j-2]+u2[i][j]) + 0.25*(u2[i][j-4]-4*u2[i][j-2]+3*u2[i][j])*(u2[i][j-4]-4*u2[i][j-2]+3*u2[i][j]);
 
 SI2 = (1.08333)*(u2[i][j-2]-2*u2[i][j]+u2[i][j+2])*(u2[i][j-2]-2*u2[i][j]+u2[i][j+2]) + 0.25*(u2[i][j-2]-u2[i][j+2])*(u2[i][j-2]-u2[i][j+2]);
 
 SI3 = (1.08333)*(u2[i][j]-2*u2[i][j+2]+u2[i][j+4])*(u2[i][j]-2*u2[i][j+2]+u2[i][j+4]) + 0.25*(3*u2[i][j]-4*u2[i][j+2]+u2[i][j+4])*(3*u2[i][j]-4*u2[i][j+2]+u2[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u2[i][j-4]- cr2*u2[i][j-2] + cr3*u2[i][j]; 
 	
 S2 = - cr5*u2[i][j-2] + cr4*u2[i][j] + cr1*u2[i][j+2]; 
 
 S3 = - cr5*u2[i][j+4] + cr4*u2[i][j+2] + cr1*u2[i][j]; 
 	
 	
 S_3 = cr1*u2[i][j+4]- cr2*u2[i][j+2] + cr3*u2[i][j];
 
 S_2 = - cr5*u2[i][j+2] + cr4*u2[i][j] + cr1*u2[i][j-2]; 
 
 S_1 = - cr5*u2[i][j-4] + cr4*u2[i][j-2] + cr1*u2[i][j]; 
 
 
                                                           
  velocity_yL[i][j+1] =S1*w1+S2*w2+S3*w3;	
 	
  velocity_yR[i][j-1] =S_1*w_1+S_2*w_2+S_3*w_3;	
 	


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u2[i-4][j]-2*u2[i-2][j]+u2[i][j])*(u2[i-4][j]-2*u2[i-2][j]+u2[i][j]) + 0.25*(u2[i-4][j]-4*u2[i-2][j]+3*u2[i][j])*(u2[i-4][j]-4*u2[i-2][j]+3*u2[i][j]);
 
 SI2 = (1.08)*(u2[i-2][j]-2*u2[i][j]+u2[i+2][j])*(u2[i-2][j]-2*u2[i][j]+u2[i+2][j]) + 0.25*(u2[i-2][j]-u2[i+2][j])*(u2[i-2][j]-u2[i+2][j]);
 
 SI3 = (1.08)*(u2[i][j]-2*u2[i+2][j]+u2[i+4][j])*(u2[i][j]-2*u2[i+2][j]+u2[i+4][j]) + 0.25*(3*u2[i][j]-4*u2[i+2][j]+u2[i+4][j])*(3*u2[i][j]-4*u2[i+2][j]+u2[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u2[i-4][j]- cr2*u2[i-2][j] + cr3*u2[i][j]; 
 	
 S2 = - cr5*u2[i-2][j] + cr4*u2[i][j] + cr1*u2[i+2][j]; 
 
 S3 = - cr5*u2[i+4][j] + cr4*u2[i+2][j] + cr1*u2[i][j]; 
 	
 	
 S_3 = cr1*u2[i+4][j]- cr2*u2[i+2][j] + cr3*u2[i][j];
 
 S_2 = - cr5*u2[i+2][j] + cr4*u2[i][j] + cr1*u2[i-2][j]; 
 
 S_1 = - cr5*u2[i-4][j] + cr4*u2[i-2][j] + cr1*u2[i][j]; 
 


v_yL[i+1][j] =S1*w1+S2*w2+S3*w3;	
 	
  v_yR[i-1][j] =S_1*w_1+S_2*w_2+S_3*w_3;	
 
 
 	 }
 	// return;
}
}


void WENO4(double u3[][points])
{

for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 
 SI1 = (1.0833)*(u3[i][j-4]-2*u3[i][j-2]+u3[i][j])*(u3[i][j-4]-2*u3[i][j-2]+u3[i][j]) + 0.25*(u3[i][j-4]-4*u3[i][j-2]+3*u3[i][j])*(u3[i][j-4]-4*u3[i][j-2]+3*u3[i][j]);
 
 SI2 = (1.08333)*(u3[i][j-2]-2*u3[i][j]+u3[i][j+2])*(u3[i][j-2]-2*u3[i][j]+u3[i][j+2]) + 0.25*(u3[i][j-2]-u3[i][j+2])*(u3[i][j-2]-u3[i][j+2]);
 
 SI3 = (1.08333)*(u3[i][j]-2*u3[i][j+2]+u3[i][j+4])*(u3[i][j]-2*u3[i][j+2]+u3[i][j+4]) + 0.25*(3*u3[i][j]-4*u3[i][j+2]+u3[i][j+4])*(3*u3[i][j]-4*u3[i][j+2]+u3[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u3[i][j-4]- cr2*u3[i][j-2] + cr3*u3[i][j]; 
 	
 S2 = - cr5*u3[i][j-2] + cr4*u3[i][j] + cr1*u3[i][j+2]; 
 
 S3 = - cr5*u3[i][j+4] + cr4*u3[i][j+2] + cr1*u3[i][j]; 
 	
 	
 S_3 = cr1*u3[i][j+4]- cr2*u3[i][j+2] + cr3*u3[i][j];
 
 S_2 = - cr5*u3[i][j+2] + cr4*u3[i][j] + cr1*u3[i][j-2]; 
 
 S_1 = - cr5*u3[i][j-4] + cr4*u3[i][j-2] + cr1*u3[i][j]; 
 
 

                                                           
  ethalpy_L[i][j+1] = S1*w1+S2*w2+S3*w3;	
 	
  ethalpy_R[i][j-1] = S_1*w_1+S_2*w_2+S_3*w_3;	
 	
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u3[i-4][j]-2*u3[i-2][j]+u3[i][j])*(u3[i-4][j]-2*u3[i-2][j]+u3[i][j]) + 0.25*(u3[i-4][j]-4*u3[i-2][j]+3*u3[i][j])*(u3[i-4][j]-4*u3[i-2][j]+3*u3[i][j]);
 
 SI2 = (1.08)*(u3[i-2][j]-2*u3[i][j]+u3[i+2][j])*(u3[i-2][j]-2*u3[i][j]+u3[i+2][j]) + 0.25*(u3[i-2][j]-u3[i+2][j])*(u3[i-2][j]-u3[i+2][j]);
 
 SI3 = (1.08)*(u3[i][j]-2*u3[i+2][j]+u3[i+4][j])*(u3[i][j]-2*u3[i+2][j]+u3[i+4][j]) + 0.25*(3*u3[i][j]-4*u3[i+2][j]+u3[i+4][j])*(3*u3[i][j]-4*u3[i+2][j]+u3[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u3[i-4][j]- cr2*u3[i-2][j] + cr3*u3[i][j]; 
 	
 S2 = - cr5*u3[i-2][j] + cr4*u3[i][j] + cr1*u3[i+2][j]; 
 
 S3 = - cr5*u3[i+4][j] + cr4*u3[i+2][j] + cr1*u3[i][j]; 
 	
 	
 S_3 = cr1*u3[i+4][j]- cr2*u3[i+2][j] + cr3*u3[i][j];
 
 S_2 = - cr5*u3[i+2][j] + cr4*u3[i][j] + cr1*u3[i-2][j]; 
 
 S_1 = - cr5*u3[i-4][j] + cr4*u3[i-2][j] + cr1*u3[i][j]; 
 
	
e_L[i+1][j] = S1*w1+S2*w2+S3*w3;	
 	
  e_R[i-1][j] = S_1*w_1+S_2*w_2+S_3*w_3;	
 
 	 
 	 }
 //return;
}
}
void WENO5(double u4[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = (1.0833)*(u4[i][j-4]-2*u4[i][j-2]+u4[i][j])*(u4[i][j-4]-2*u4[i][j-2]+u4[i][j]) + 0.25*(u4[i][j-4]-4*u4[i][j-2]+3*u4[i][j])*(u4[i][j-4]-4*u4[i][j-2]+3*u4[i][j]);
 
 SI2 = (1.08333)*(u4[i][j-2]-2*u4[i][j]+u4[i][j+2])*(u4[i][j-2]-2*u4[i][j]+u4[i][j+2]) + 0.25*(u4[i][j-2]-u4[i][j+2])*(u4[i][j-2]-u4[i][j+2]);
 
 SI3 = (1.08333)*(u4[i][j]-2*u4[i][j+2]+u4[i][j+4])*(u4[i][j]-2*u4[i][j+2]+u4[i][j+4]) + 0.25*(3*u4[i][j]-4*u4[i][j+2]+u4[i][j+4])*(3*u4[i][j]-4*u4[i][j+2]+u4[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u4[i][j-4]- cr2*u4[i][j-2] + cr3*u4[i][j]; 
 	
 S2 = - cr5*u4[i][j-2] + cr4*u4[i][j] + cr1*u4[i][j+2]; 
 
 S3 = - cr5*u4[i][j+4] + cr4*u4[i][j+2] + cr1*u4[i][j]; 
 	
 	
 S_3 = cr1*u4[i][j+4]- cr2*u4[i][j+2] + cr3*u4[i][j];
 
 S_2 = - cr5*u4[i][j+2] + cr4*u4[i][j] + cr1*u4[i][j-2]; 
 
 S_1 = - cr5*u4[i][j-4] + cr4*u4[i][j-2] + cr1*u4[i][j]; 
 
  
 
  velocity_xL[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
  velocity_xR[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u4[i-4][j]-2*u4[i-2][j]+u4[i][j])*(u4[i-4][j]-2*u4[i-2][j]+u4[i][j]) + 0.25*(u4[i-4][j]-4*u4[i-2][j]+3*u4[i][j])*(u4[i-4][j]-4*u4[i-2][j]+3*u4[i][j]);
 
 SI2 = (1.08)*(u4[i-2][j]-2*u4[i][j]+u4[i+2][j])*(u4[i-2][j]-2*u4[i][j]+u4[i+2][j]) + 0.25*(u4[i-2][j]-u4[i+2][j])*(u4[i-2][j]-u4[i+2][j]);
 
 SI3 = (1.08)*(u4[i][j]-2*u4[i+2][j]+u4[i+4][j])*(u4[i][j]-2*u4[i+2][j]+u4[i+4][j]) + 0.25*(3*u4[i][j]-4*u4[i+2][j]+u4[i+4][j])*(3*u4[i][j]-4*u4[i+2][j]+u4[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u4[i-4][j]- cr2*u4[i-2][j] + cr3*u4[i][j]; 
 	
 S2 = - cr5*u4[i-2][j] + cr4*u4[i][j] + cr1*u4[i+2][j]; 
 
 S3 = - cr5*u4[i+4][j] + cr4*u4[i+2][j] + cr1*u4[i][j]; 
 	
 	
 S_3 = cr1*u4[i+4][j]- cr2*u4[i+2][j] + cr3*u4[i][j];
 
 S_2 = - cr5*u4[i+2][j] + cr4*u4[i][j] + cr1*u4[i-2][j]; 
 
 S_1 = - cr5*u4[i-4][j] + cr4*u4[i-2][j] + cr1*u4[i][j]; 
 
 
 
 v_xL[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
  v_xR[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
 	 
 	 }
//return;

}
}
 
 void WENO6(double u5[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = (1.0833)*(u5[i][j-4]-2*u5[i][j-2]+u5[i][j])*(u5[i][j-4]-2*u5[i][j-2]+u5[i][j]) + 0.25*(u5[i][j-4]-4*u5[i][j-2]+3*u5[i][j])*(u5[i][j-4]-4*u5[i][j-2]+3*u5[i][j]);
 
 SI2 = (1.08333)*(u5[i][j-2]-2*u5[i][j]+u5[i][j+2])*(u5[i][j-2]-2*u5[i][j]+u5[i][j+2]) + 0.25*(u5[i][j-2]-u5[i][j+2])*(u5[i][j-2]-u5[i][j+2]);
 
 SI3 = (1.08333)*(u5[i][j]-2*u5[i][j+2]+u5[i][j+4])*(u5[i][j]-2*u5[i][j+2]+u5[i][j+4]) + 0.25*(3*u5[i][j]-4*u5[i][j+2]+u5[i][j+4])*(3*u5[i][j]-4*u5[i][j+2]+u5[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u5[i][j-4]- cr2*u5[i][j-2] + cr3*u5[i][j]; 
 	
 S2 = - cr5*u5[i][j-2] + cr4*u5[i][j] + cr1*u5[i][j+2]; 
 
 S3 = - cr5*u5[i][j+4] + cr4*u5[i][j+2] + cr1*u5[i][j]; 
 	
 	
 S_3 = cr1*u5[i][j+4]- cr2*u5[i][j+2] + cr3*u5[i][j];
 
 S_2 = - cr5*u5[i][j+2] + cr4*u5[i][j] + cr1*u5[i][j-2]; 
 
 S_1 = - cr5*u5[i][j-4] + cr4*u5[i][j-2] + cr1*u5[i][j]; 

 
  magnetic_L[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
  magnetic_R[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u5[i-4][j]-2*u5[i-2][j]+u5[i][j])*(u5[i-4][j]-2*u5[i-2][j]+u5[i][j]) + 0.25*(u5[i-4][j]-4*u5[i-2][j]+3*u5[i][j])*(u5[i-4][j]-4*u5[i-2][j]+3*u5[i][j]);
 
 SI2 = (1.08)*(u5[i-2][j]-2*u5[i][j]+u5[i+2][j])*(u5[i-2][j]-2*u5[i][j]+u5[i+2][j]) + 0.25*(u5[i-2][j]-u5[i+2][j])*(u5[i-2][j]-u5[i+2][j]);
 
 SI3 = (1.08)*(u5[i][j]-2*u5[i+2][j]+u5[i+4][j])*(u5[i][j]-2*u5[i+2][j]+u5[i+4][j]) + 0.25*(3*u5[i][j]-4*u5[i+2][j]+u5[i+4][j])*(3*u5[i][j]-4*u5[i+2][j]+u5[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u5[i-4][j]- cr2*u5[i-2][j] + cr3*u5[i][j]; 
 	
 S2 = - cr5*u5[i-2][j] + cr4*u5[i][j] + cr1*u5[i+2][j]; 
 
 S3 = - cr5*u5[i+4][j] + cr4*u5[i+2][j] + cr1*u5[i][j]; 
 	
 	
 S_3 = cr1*u5[i+4][j]- cr2*u5[i+2][j] + cr3*u5[i][j];
 
 S_2 = - cr5*u5[i+2][j] + cr4*u5[i][j] + cr1*u5[i-2][j]; 
 
 S_1 = - cr5*u5[i-4][j] + cr4*u5[i-2][j] + cr1*u5[i][j]; 
 
 
  BY_L[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
 BY_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
 	 
 	 }


}
}
void WENO7(double u6[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = (1.0833)*(u6[i][j-4]-2*u6[i][j-2]+u6[i][j])*(u6[i][j-4]-2*u6[i][j-2]+u6[i][j]) + 0.25*(u6[i][j-4]-4*u6[i][j-2]+3*u6[i][j])*(u6[i][j-4]-4*u6[i][j-2]+3*u6[i][j]);
 
 SI2 = (1.08333)*(u6[i][j-2]-2*u6[i][j]+u6[i][j+2])*(u6[i][j-2]-2*u6[i][j]+u6[i][j+2]) + 0.25*(u6[i][j-2]-u6[i][j+2])*(u6[i][j-2]-u6[i][j+2]);
 
 SI3 = (1.08333)*(u6[i][j]-2*u6[i][j+2]+u6[i][j+4])*(u6[i][j]-2*u6[i][j+2]+u6[i][j+4]) + 0.25*(3*u6[i][j]-4*u6[i][j+2]+u6[i][j+4])*(3*u6[i][j]-4*u6[i][j+2]+u6[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u6[i][j-4]- cr2*u6[i][j-2] + cr3*u6[i][j]; 
 	
 S2 = - cr5*u6[i][j-2] + cr4*u6[i][j] + cr1*u6[i][j+2]; 
 
 S3 = - cr5*u6[i][j+4] + cr4*u6[i][j+2] + cr1*u6[i][j]; 
 	
 	
 S_3 = cr1*u6[i][j+4]- cr2*u6[i][j+2] + cr3*u6[i][j];
 
 S_2 = - cr5*u6[i][j+2] + cr4*u6[i][j] + cr1*u6[i][j-2]; 
 
 S_1 = - cr5*u6[i][j-4] + cr4*u6[i][j-2] + cr1*u6[i][j]; 

 
  bx_L[i][j+1] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
  bx_R[i][j-1] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u6[i-4][j]-2*u6[i-2][j]+u6[i][j])*(u6[i-4][j]-2*u6[i-2][j]+u6[i][j]) + 0.25*(u6[i-4][j]-4*u6[i-2][j]+3*u6[i][j])*(u6[i-4][j]-4*u6[i-2][j]+3*u6[i][j]);
 
 SI2 = (1.08)*(u6[i-2][j]-2*u6[i][j]+u6[i+2][j])*(u6[i-2][j]-2*u6[i][j]+u6[i+2][j]) + 0.25*(u6[i-2][j]-u6[i+2][j])*(u6[i-2][j]-u6[i+2][j]);
 
 SI3 = (1.08)*(u6[i][j]-2*u6[i+2][j]+u6[i+4][j])*(u6[i][j]-2*u6[i+2][j]+u6[i+4][j]) + 0.25*(3*u6[i][j]-4*u6[i+2][j]+u6[i+4][j])*(3*u6[i][j]-4*u6[i+2][j]+u6[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u6[i-4][j]- cr2*u6[i-2][j] + cr3*u6[i][j]; 
 	
 S2 = - cr5*u6[i-2][j] + cr4*u6[i][j] + cr1*u6[i+2][j]; 
 
 S3 = - cr5*u6[i+4][j] + cr4*u6[i+2][j] + cr1*u6[i][j]; 
 	
 	
 S_3 = cr1*u6[i+4][j]- cr2*u6[i+2][j] + cr3*u6[i][j];
 
 S_2 = - cr5*u6[i+2][j] + cr4*u6[i][j] + cr1*u6[i-2][j]; 
 
 S_1 = - cr5*u6[i-4][j] + cr4*u6[i-2][j] + cr1*u6[i][j]; 
 
 
  BX_L[i+1][j] = (S1*w1)+(S2*w2)+(S3*w3);	
 	
 BX_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
 	 }

}
return;
}

void WENO8(double u7[][points])
{

for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
{
    
 // if(i==10)
  //cin>>j;
 
 SI1 = (1.0833)*(u7[i][j-4]-2*u7[i][j-2]+u7[i][j])*(u7[i][j-4]-2*u7[i][j-2]+u7[i][j]) + 0.25*(u7[i][j-4]-4*u7[i][j-2]+3*u7[i][j])*(u7[i][j-4]-4*u7[i][j-2]+3*u7[i][j]);
 
 SI2 = (1.08333)*(u7[i][j-2]-2*u7[i][j]+u7[i][j+2])*(u7[i][j-2]-2*u7[i][j]+u7[i][j+2]) + 0.25*(u7[i][j-2]-u7[i][j+2])*(u7[i][j-2]-u7[i][j+2]);
 
 SI3 = (1.08333)*(u7[i][j]-2*u7[i][j+2]+u7[i][j+4])*(u7[i][j]-2*u7[i][j+2]+u7[i][j+4]) + 0.25*(3*u7[i][j]-4*u7[i][j+2]+u7[i][j+4])*(3*u7[i][j]-4*u7[i][j+2]+u7[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u7[i][j-4]- cr2*u7[i][j-2] + cr3*u7[i][j]; 
 	
 S2 = - cr5*u7[i][j-2] + cr4*u7[i][j] + cr1*u7[i][j+2]; 
 
 S3 = - cr5*u7[i][j+4] + cr4*u7[i][j+2] + cr1*u7[i][j]; 
 	
 	
 S_3 = cr1*u7[i][j+4]- cr2*u7[i][j+2] + cr3*u7[i][j];
 
 S_2 = - cr5*u7[i][j+2] + cr4*u7[i][j] + cr1*u7[i][j-2]; 
 
 S_1 = - cr5*u7[i][j-4] + cr4*u7[i][j-2] + cr1*u7[i][j]; 
 
 
                                                           
  velocity_zL[i][j+1] = S1*w1+S2*w2+S3*w3;	
 	
  velocity_zR[i][j-1] = S_1*w_1+S_2*w_2+S_3*w_3;	
 	


 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u7[i-4][j]-2*u7[i-2][j]+u7[i][j])*(u7[i-4][j]-2*u7[i-2][j]+u7[i][j]) + 0.25*(u7[i-4][j]-4*u7[i-2][j]+3*u7[i][j])*(u7[i-4][j]-4*u7[i-2][j]+3*u7[i][j]);
 
 SI2 = (1.08)*(u7[i-2][j]-2*u7[i][j]+u7[i+2][j])*(u7[i-2][j]-2*u7[i][j]+u7[i+2][j]) + 0.25*(u7[i-2][j]-u7[i+2][j])*(u7[i-2][j]-u7[i+2][j]);
 
 SI3 = (1.08)*(u7[i][j]-2*u7[i+2][j]+u7[i+4][j])*(u7[i][j]-2*u7[i+2][j]+u7[i+4][j]) + 0.25*(3*u7[i][j]-4*u7[i+2][j]+u7[i+4][j])*(3*u7[i][j]-4*u7[i+2][j]+u7[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u7[i-4][j]- cr2*u7[i-2][j] + cr3*u7[i][j]; 
 	
 S2 = - cr5*u7[i-2][j] + cr4*u7[i][j] + cr1*u7[i+2][j]; 
 
 S3 = - cr5*u7[i+4][j] + cr4*u7[i+2][j] + cr1*u7[i][j]; 
 	
 	
 S_3 = cr1*u7[i+4][j]- cr2*u7[i+2][j] + cr3*u7[i][j];
 
 S_2 = - cr5*u7[i+2][j] + cr4*u7[i][j] + cr1*u7[i-2][j]; 
 
 S_1 = - cr5*u7[i-4][j] + cr4*u7[i-2][j] + cr1*u7[i][j]; 
 


v_zL[i+1][j] = S1*w1+S2*w2+S3*w3;	
 	
  v_zR[i-1][j] = S_1*w_1+S_2*w_2+S_3*w_3;	
 
 
 	 }
 	// return;
}
}

void WENO9(double u8[][points]) //,double u1[],double u2[], double u3[])
 {
 	
for(i=4; i<=size-4; i=i+2)
 {
 	
for(j=4; j<=size-4; j=j+2)
 {
  
 // if(i==10)
  //cin>>j;
 
SI1 = (1.0833)*(u8[i][j-4]-2*u8[i][j-2]+u8[i][j])*(u8[i][j-4]-2*u8[i][j-2]+u8[i][j]) + 0.25*(u8[i][j-4]-4*u8[i][j-2]+3*u8[i][j])*(u8[i][j-4]-4*u8[i][j-2]+3*u8[i][j]);
 
 SI2 = (1.08333)*(u8[i][j-2]-2*u8[i][j]+u8[i][j+2])*(u8[i][j-2]-2*u8[i][j]+u8[i][j+2]) + 0.25*(u8[i][j-2]-u8[i][j+2])*(u8[i][j-2]-u8[i][j+2]);
 
 SI3 = (1.08333)*(u8[i][j]-2*u8[i][j+2]+u8[i][j+4])*(u8[i][j]-2*u8[i][j+2]+u8[i][j+4]) + 0.25*(3*u8[i][j]-4*u8[i][j+2]+u8[i][j+4])*(3*u8[i][j]-4*u8[i][j+2]+u8[i][j+4]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u8[i][j-4]- cr2*u8[i][j-2] + cr3*u8[i][j]; 
 	
 S2 = - cr5*u8[i][j-2] + cr4*u8[i][j] + cr1*u8[i][j+2]; 
 
 S3 = - cr5*u8[i][j+4] + cr4*u8[i][j+2] + cr1*u8[i][j]; 
 	
 	
 S_3 = cr1*u8[i][j+4]- cr2*u8[i][j+2] + cr3*u8[i][j];
 
 S_2 = - cr5*u8[i][j+2] + cr4*u8[i][j] + cr1*u8[i][j-2]; 
 
 S_1 = - cr5*u8[i][j-4] + cr4*u8[i][j-2] + cr1*u8[i][j]; 

 
  bz_L[i][j+1] =  (S1*w1)+(S2*w2)+(S3*w3);	
 	
  bz_R[i][j-1] =  (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
    
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////// for y-direction	
 	
  SI1 = (1.08)*(u8[i-4][j]-2*u8[i-2][j]+u8[i][j])*(u8[i-4][j]-2*u8[i-2][j]+u8[i][j]) + 0.25*(u8[i-4][j]-4*u8[i-2][j]+3*u8[i][j])*(u8[i-4][j]-4*u8[i-2][j]+3*u8[i][j]);
 
 SI2 = (1.08)*(u8[i-2][j]-2*u8[i][j]+u8[i+2][j])*(u8[i-2][j]-2*u8[i][j]+u8[i+2][j]) + 0.25*(u8[i-2][j]-u8[i+2][j])*(u8[i-2][j]-u8[i+2][j]);
 
 SI3 = (1.08)*(u8[i][j]-2*u8[i+2][j]+u8[i+4][j])*(u8[i][j]-2*u8[i+2][j]+u8[i+4][j]) + 0.25*(3*u8[i][j]-4*u8[i+2][j]+u8[i+4][j])*(3*u8[i][j]-4*u8[i+2][j]+u8[i+4][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u8[i-4][j]- cr2*u8[i-2][j] + cr3*u8[i][j]; 
 	
 S2 = - cr5*u8[i-2][j] + cr4*u8[i][j] + cr1*u8[i+2][j]; 
 
 S3 = - cr5*u8[i+4][j] + cr4*u8[i+2][j] + cr1*u8[i][j]; 
 	
 	
 S_3 = cr1*u8[i+4][j]- cr2*u8[i+2][j] + cr3*u8[i][j];
 
 S_2 = - cr5*u8[i+2][j] + cr4*u8[i][j] + cr1*u8[i-2][j]; 
 
 S_1 = - cr5*u8[i-4][j] + cr4*u8[i-2][j] + cr1*u8[i][j]; 
 
 
  BZ_L[i+1][j] =  (S1*w1)+(S2*w2)+(S3*w3);	
 	
 BZ_R[i-1][j] = (S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
 	
	 //cout<<BZ_L[i+1][j]<<endl;
	 
	  }

}
return;
}

