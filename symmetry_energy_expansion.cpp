//
//  NStoHIC
//
//  Created by Nanxi Yao on 2/26/22.
//


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <math.h>


using namespace std;

const double nsat =0.161;
const double hc = 197.3;
const double mN=938.0;
const double hbar=6.582119569*pow(10.0,-22.0);
const double c=2.998*pow(10.0,23.0);
const double mmu=105.7/c/c;
const double me=0.51099895/c/c;

const double PI  =3.141592653589793238463;

double getFermiEnergy(double mass, double ni)
{
    double kF=pow(ni/2.0*6*PI*PI,1.0/3.0)*hbar;
    double epsilon=pow(mass,4.0)*pow(c,5.0)/PI/PI/hbar/hbar/hbar;
    double x=kF/mass/c;
    double eden=epsilon/8.0*((2*pow(x,3.0)+x)*pow((1+x*x),1.0/2.0)-asinh(x));
    return eden;
    
}
double defYQ(double nB, double E, double L, double K, double J)
{
    double x=(nB-nsat)/3.0/nsat;
    double esym=(E+L*x+K*x*x/2.0+J*x*x*x/6.0)/hc;
    double term2=288*pow(esym,12.0)*pow(nB,2.0)+PI*PI*pow(esym,9.0)*pow(nB,3.0);
    double term1=-24*pow(esym,6.0)*nB+pow(2,1.0/2.0)*pow(term2,1.0/2.0);
    double term3=pow(term1,1.0/3.0);
    double yq=1.0/16.0*(8-pow(PI,4.0/3.0)*nB/term3/pow(2,1.0/3.0)+pow(PI/2,2.0/3.0)*term3/pow(esym,3.0));
  //  double yq=64.0/3.0/PI/PI/nsat/(3.0*x+1.0)*power*power*power;
    
    double powernsat = E/hc;
 //   double yqnsat =64.0/3.0/PI/PI/nsat*powernsat*powernsat*powernsat;
  //  if(esym<0)  return 0;
    
   // else
        return yq;
}

// This conversion only considers expanding from NS to symmetric matter
double getEHIC(double eNS,double rho,double e0,double L, double K, double J, double yQ)
{
    double delta=0;
    return eNS-(e0+L*(rho/nsat-1)/3+K/18*(rho/nsat-1)*(rho/nsat-1)+J/162*((rho/nsat-1)*(rho/nsat-1)*(rho/nsat-1)))*(1-2*yQ)*(1-2*yQ)*rho;
}

//expanding from NS to symmetric matter and then to HIC (y_HIC varies)
double getEHIC2(double eNS, double rho,double e0,double L, double K, double J, double Y_HIC, double yQ)
{
    return eNS-(e0+L*(rho/nsat-1)/3+K/18*(rho/nsat-1)*(rho/nsat-1)+J/162*((rho/nsat-1)*(rho/nsat-1)*(rho/nsat-1)))*4*((Y_HIC-yQ)+(yQ*yQ-Y_HIC*Y_HIC))*rho;
}

//checks if converted cs2 is causal for nB<rholimit, number of acausal points less than numlimit
bool isCausal(double *cs2, double rhoB[], double numlimit, int limit)
{
    int count=0;
    for(int i=0;i<limit;i++)
    {
        if (rhoB[i]/nsat>0 && rhoB[i]/nsat<0.9)
        {
            
            if(cs2[i]>1)
            {
                count++;
            }
        }
        else if(rhoB[i]/nsat>0.9 )
        {
            if(cs2[i]>1 || cs2[i]<0)
            {
                count++;
            }
        }
        
        
    }
    if(count > numlimit)
    {
        return false;
    }
    return true;
    
}


double getnumericalderiv(double x[], double y[], int i)
{
    double deriv = (y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
    return deriv;
}


double** convert(double *&rhoB, double *&energy, double e0, double L, double K, double J, int limit, double YHIC, double switch_pt)
{
    double enerHIC[limit], eden[limit], pressure[limit],cs2[limit],eoverrho[limit],rhoB2[limit];
    double *ne =new double[limit];
    double *energysub=new double[limit];
    double *yQ_arr = new double[limit];

    int start=0;
    int end=0;
    for(int i=0;i<limit;i++)
    {
     //   double yQ=YQ_cvgt(switch_pt, e0, L, K, rhoB[i]);
        double yQ=defYQ(rhoB[i], e0, L, K,J);
        yQ_arr[i]=yQ;
      
    }
    for(int i=0;i<limit;i++)
    {
        double yQ=yQ_arr[i];
        ne[i]=rhoB[i]*yQ;
        if(ne[i]<0) ne[i]=-ne[i];
        double neval=ne[i];
        energysub[i]=energy[i]-getFermiEnergy(me,neval);
        eden[i]=getEHIC2(energysub[i],rhoB[i],e0,L,K,J, YHIC,yQ);
       // eden[i]=getEHIC(energy[i],rhoB[i],e0,L,K,J,YQ[i]);
       // eden[i]=getEHIC(energy[i],rhoB[i],e0,L,K,J,0.1);
        eoverrho[i]=eden[i]/rhoB[i];
        rhoB2[i]=rhoB[i];
    }
    for(int i=0;i<limit;i++)
    {
        pressure[i]=getnumericalderiv(rhoB2,eoverrho,i)*rhoB2[i]*rhoB2[i];
    }

    for(int i=0;i<limit;i++)
    {
        cs2[i]=getnumericalderiv(eden,pressure,i);
    }
    double** rhoB_cs2 = 0;
    rhoB_cs2 = new double*[4];
    rhoB_cs2[0] = new double[limit];
    rhoB_cs2[1] = new double[limit];
    rhoB_cs2[2] = new double[limit];
    rhoB_cs2[3] = new double[limit];
    for (int w = 0; w < limit; w++)
        {
        rhoB_cs2[0][w]=rhoB[w];
    //   rhoB_cs2[1][w]=yQ_arr[w];
    //  rhoB_cs2[1][w]=fermi[w];
     　  rhoB_cs2[3][w]=cs2[w];
   //   cout<<w<<" "<<rhoB2[w]<<endl;
      rhoB_cs2[2][w]=pressure[w];
      rhoB_cs2[1][w]=eden[w];
    //  rhoB_cs2[1][w]=eoverrho[w];
                }
  /*  for(int i=0;i<limit;i++)
    {
    //    cout<<i<<" "<<cs2<<endl;
    }*/
    return rhoB_cs2;
  //  return cs2;
}


int main() {
 
    FILE *finput = fopen("eos1.txt","r");
    int limit=6075;             //**Open the input file.
    double *energy = new double[limit];
    double *muB = new double[limit];
    double *rhoB = new double[limit];
    double *pressure =new double[limit];
    double drop2;

    for(int i=0;i<limit;i++)
    {
        fscanf(finput,"%lf\t%lf\t%lf\t%lf\t%lf\n",&rhoB[i],&pressure[i],&energy[i],&drop2,&muB[i]);
    }
    fclose(finput);
    /*get original EOS info*/
 /*   ofstream myfile2("original_eos4_cs2.txt");
   //  ofstream myfile2("sly4_allyq.txt");
  //
    double oldcs2[limit];
    double allyQ[limit];
    for(int i=0;i<limit;i++)
    {
          oldcs2[i]=getnumericalderiv(energy,pressure,i);
       //  myfile<<rhoB[i]/nsat<<" "<<energy[i]/(rhoB[i]*nsat)<<endl;
      //  myfile2<<rhoB[i]/nsat<<" "<<oldcs2[i]<<endl;
    }*/
    
     
    //get yQ information
    
   /*  double e0=31.7;
       double L=58.7;
       double K=-100;
       double J=32;
    ofstream myfile2("getYQ.txt");
    for(int i=0;i<limit;i++)
    {
        myfile2<<rhoB[i]/nsat<<" "<<getgausYQ(e0, L, K,rhoB[i])<<endl;
    }
*/
    
    /*Loop through different coefficients*/
  
  　 int count=0;
    double e0=27.5;
    double L=30;
    double K=-220;
    double J=-200;
    vector<vector<double>> all_cs2;
    vector<vector<double>> all_rhoB;
    vector<vector<double>> all_eden;
    vector<vector<double>> all_pressure;
    ofstream myfile("all_cs2_eden_pressure_eos1.txt");
   // ofstream myfile("all_yQ.txt");
   // ofstream myfile2("coefficients_eos3_nocausallim.txt");
    int i1lim=13;
    int i2lim=10;
    int i3lim=8;
    int i4lim=9;
    bool causal=true;
    for(int i1=0;i1<i1lim;i1++)
    {
        for(int i2=0;i2<i2lim;i2++)
        {
           for(int i3=0;i3<i3lim;i3++)
           {
              for(int i4=0;i4<i4lim;i4++)
                {
                    double x=0;
                    double esym=0;
                    bool reject=false;
                    for(int i5=0;i5<limit;i5++)
                    {
                    x=(rhoB[i5]-nsat)/3.0/nsat;
                    esym=(e0+L*x+K*x*x/2.0+J*x*x*x/6.0)/hc;
                    if(esym<0 && rhoB[i5] < 6*nsat)
                        {
                            reject=true;
                            break;
                        }
                    }
                 
                    if(esym<0 && reject == true)
                    {
                        J=J+100;
                        continue;
                    }
                  
                    else
                    {
                    double **cs2 = convert(rhoB,energy,e0,L,K,J,limit,0.5,0.15);
                    bool causal=isCausal(cs2[3],rhoB,10,limit);
                   // cout<<causal<<endl;
                    vector<double> single_cs2;
                    vector<double> single_rhoB;
                    vector<double> single_eden;
                    vector<double> single_pressure;
                   if(causal)
                    {
                        for(int k=0;k<limit;k++)
                        {
                            single_eden.push_back(cs2[1][k]);
                            single_rhoB.push_back(cs2[0][k]);
                            single_pressure.push_back(cs2[2][k]);
                            single_cs2.push_back(cs2[3][k]);
                        //    cout<<cs2[k]<<" ";
                        }
                        all_cs2.push_back(single_cs2);
                        all_rhoB.push_back(single_rhoB);
                        all_eden.push_back(single_eden);
                        all_pressure.push_back(single_pressure);
                      //  myfile2<<e0<<" "<<L<<" "<<K<<" "<<J<<endl;
                       
                       
                        count++;
                    }
                    }
           
                J=J+100;
              
                 
                }
                J=-200;
                K=K+50;
          
            }
            K=-220;
            L=L+10;
        }
        L=30;
        e0=e0+1;
        cout<<e0<<" "<<count<<endl;
    }
    //output cs2 information
   for(int i=0;i<count;i++)
    {
        for(int j=0;j<limit;j++)
        {
          
            myfile<<all_rhoB[i][j]<<" "<<all_eden[i][j]<<" "<<all_pressure[i][j]<<" "<<all_cs2[i][j]<<endl;
       
        }
    }

//get single conversino information with specific symmetry energy coefficients
  /*    double e0=32.0;
      double L=47.46;
      double K=-115.13;
      double J=0;*/


    
 double **cs2=convert(rhoB,energy,e0,L,K,J,limit,0.5,0.15);
    ofstream myfile("eos1_4nsat.txt");
    int countt=0;
 //   ofstream myfile("nqmuq_eos3_01.txt");
    for(int i=0;i<limit;i++)
    {
      
            myfile<<cs2[0][i]/nsat<<" "<<cs2[1][i]<<endl;
      //  myfile<<rhoB[i]/nsat<<" "<<YQ_cvgt(0.15, e0, L, K, rhoB[i])<<endl;
      //  myfile<<rhoB[i]/nsat<<" "<<YQ_cvgt_smooth(rhoB, 0.15, e0, L, K, limit,i)<<endl;
    }*/

    
}

