#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <istream>
#include <cstring>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <stdlib.h>
using namespace std;

typedef complex<double> comp;
#ifndef M_PI
#define M_PI    3.1415926535897
#endif
//-----PHYSICAL CONSTANTS (QCD-SYS->c=1,hbar=1,4pi*e_0=1)---------//
const double m=4.14/2;                                        //MASS for evolution
const double hbar=1;                                          // h/2pi
const double T = 0.4;                                         //Temperature
const short int N_c = 3;                                      //# of colors
const short int N_f = 2;                                      //# quark flavours
const double C_f = (double)(N_c*N_c-1)/((double)(2*N_c));
const double alpha_s = 0.41;							      //Strong coupling
double alpha_s_im = 0.41;                                     //this value will run in auto procedure
const double gamma_E = 0.57721566490153286060;                //euler constant
double g = sqrt(4*M_PI*alpha_s_im);
double m_D = pow(g*T,2)/3. *(N_c+(double)N_f/2.);			  //Debye Mass
//----------------Grid--------------------------
const size_t N=2000;                                           //Gridpoints
double A = 0.05;                                               //Spatial Step
comp EPS = 0.0001;                                             //Time Step
int TIME =250000;                                              //Total Time
int t=0;
double L = (static_cast <double> (N))*A;                       //Box Size
comp alpha = EPS/(4.*A*A*m);
//----------------------------------------------
char buf[255],initbuf[255];
comp normalization,energy,energy1,energy0;
int POTENTIAL,INITCONDITION,SYMM,RECORD,AUTO;
clock_t update_time;
int record_point=2,snap_num=4;
const int snaptime = TIME/(snap_num);
comp im = comp(0,1);
//-----------FILES-------------------------
ofstream debug ("data/INFO.txt");
ofstream pot("data/pot.txt");
ofstream ener1("data/energies1.txt"); //files for qq-energies
//----------Arrays and Matrices-------------
vector <comp> diag(N);                                        //holds diagnoal
vector <comp> b(N);                                           //update vactor
vector <comp> psi(N);                                         //holds solution
vector <comp> V(N);                                           //holds potential
vector <comp> off_diag_star(N);                               //necessary for LU_algorithm
vector <comp> b_star(N);                                      //----------"--------------
vector < vector<comp> > psi_snap(N,vector<comp>(snap_num));   //holds snapshots
cvec PSI_1(N);                                                //holds excited state

//------------Functions--------------------


//Here the Inital conditions and potential are set

void set_init(){
    //-------set initial wave---------------
double sgn,x;
switch (INITCONDITION)
{
case 1: //random
    sprintf(initbuf,"Random");
    for (int i=0;i<N;i++)
    {
        sgn=(double)1./(rand()%100 +1);
        if(sgn>0.5)
        {
           psi[i]=(double)-1./(rand()%100 +1);}
        else{
            psi[i]=(double)1./(rand()%100 +1);
        }
    }
    break;
case 2: //delta peak
    sprintf(initbuf,"Delta peak");
    for (int i=0;i<N;i++)
    {
           if(i==N/2)
           {
                 psi[N/2]=1.;
           }
           else
           {
                 psi[i]=0.0;
           }
    }
    break;
case 3://radial of coulomb
    sprintf(initbuf,"Radial coulomb");
    for (int i=0;i<N;i++)
    {
           x=i*A;
           psi[i]=2*x*exp(-x);
    }
    break;
case 4: //Particle in a box ground state
    sprintf(initbuf,"Box ground state");
    for (int i=0;i<N;i++)
    {
           x=i*A;
           psi[i]= sqrt(2./(N*A)) * sin(M_PI * x /(N*A));
    }
    break;
case 5: //advanced dirac delta
    sprintf(initbuf,"Dirac delta");
    for (int i=0;i<N;i++)
    {
           x=i*A;
           psi[i]= -6.*N_c*(2./(M_PI*A*A)) * i/(4*i*i-1) *pow(-1,i+1);
    }
    break;
default:
    break;
}

char name[255];
sprintf(name,"initial");
wfnc_out(psi,name);//save initial wave to datafile

}

void set_pot(){
    double x;
    double constant = 0.0;
  //set potential
  switch (POTENTIAL)
  {
  case 1: //coulomb
         SYMM = 1;
         sprintf(buf,"coulomb");
         for (int i=1;i<N;i++)
         {
             x=A*(double)i;
             V[i]=-pow(x,-1);
             V[i]+=(comp)constant;
         }
         V[0]=(double)-1./(A*A);
         break;
  case 2://oscillator
         SYMM = 1;
         sprintf(buf,"osc");
         for (int i=0;i<N;i++)
         {
             x = ((double)i-(double)N/2)*A;
                V[i]=0.5*x*x;
         }
         break;
  case 3://particle in a box
         SYMM = 2;
         sprintf(buf,"box");
         for (int i=0;i<N;i++)
         {
                V[i]=0.;
         }
         break;
case 4: //QQ
         SYMM = 1;
         sprintf(buf,"QQ");
         for (int i=0;i<N;i++)
         {
             x=i*A;
             V[i]=(comp)-C_f*alpha_s*pow(x,-1) +(comp)im*((C_f*(alpha_s_im*pow(x,2)*T*pow(m_D,2)))/((comp)6.) * (2*gamma_E+ log(pow(x*m_D,2))-(comp)8./3.)) + constant;
         }
         V[0]=-1./(pow(A,2));
         break;
  default:
         break;
  }

  //set diagonal

  for (int i=0;i<N;i++)
  {
         diag[i]= 1.0 + EPS/(2.0) * (1/(A*A*m) + V[i]);
  }

}

void set_param()
{
set_pot();
set_init();
}

void record_snap(vector <comp>& wfnc, int num)//saves snapshots in /data/snapshots/
{
       char name [255];
       ofstream out;

       sprintf(name,"data/snapshots/wave_%d.txt",num);
       out.open(name, ios::out);
       out.precision(10);
       switch (SYMM)
       {
       case 1:
              for (int sx =0; sx<N;sx++)
       {
              out << sx << "       " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
       }
              break;
       case 2:
              for (int sx =0; sx<N;sx++)
       {
              out << -N/2. + sx << "     " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
       }
              break;
       default:
              break;
       }

       out.close();
}

void wfnc_out(vector <comp>& wfnc, char* label) //save wavefunction in /data
{
    ofstream out;
    sprintf(label,"data/wave_%s.txt",label);
    out.open(label, ios::out);
    out.precision(10);
    switch (SYMM)
    {
    case 1:
           for (int sx =0; sx<N;sx++)
    {
           out << sx << "       " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
    }
           break;
    case 2:
           for (int sx =0; sx<N;sx++)
    {
           out << -N/2. + sx << "     " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
    }
           break;
    default:
           break;
    }

    out.close();
}
//Calculate norm of wave
inline void normalize(vector <comp>& wfnc)
{
       comp sum;
       //via trapez rule
       vector <comp> tmp (2);
       for (int i=0;i<N-1;i++)
       {
              tmp[0] = conj(wfnc[i])*wfnc[i];
              tmp[1] = conj(wfnc[i+1])*wfnc[i+1];
              sum += A*(tmp[0]+0.5*(tmp[1]-tmp[0]));
       }
       sum += -0.5*A*conj(wfnc[N-1])*wfnc[N-1];
    sum=sqrt(sum);
       for (int i=0;i<N;i++)
       {
              wfnc[i]/=sum;
       }
       normalization=sum;

       //Gives error message if norm = inf or NaN
       if (sum!=sum || sum==(comp)0)
       {
              cout << "Problem with normalization!";
       }
       tmp.clear();
}
//calculate energy
comp wfncEnergy(vector <comp>& wfnc)
{

       complex <double> res=0;
       //via trapez
       vector <comp> tmp(2);
       for (int sx=1;sx<N-2;sx++)
              {
                     tmp[0] = V[sx]*conj(wfnc[sx])*wfnc[sx]- conj(wfnc[sx])*hbar*hbar*(wfnc[sx+1] + wfnc[sx-1] - 2.*wfnc[sx])/(2.*A*A*m);
                     tmp[1] = V[sx+1]*conj(wfnc[sx+1])*wfnc[sx+1]- conj(wfnc[sx+1])*hbar*hbar*(wfnc[sx+1+1] + wfnc[sx-1+1] - 2.*wfnc[sx+1])/(2.*A*A*m);
                  res += A*(tmp[0]+0.5*(tmp[1]-tmp[0]));
              }
       res+= -0.5*A*V[N-2]*conj(wfnc[N-2])*wfnc[N-2]- conj(wfnc[N-2])*hbar*hbar*(wfnc[N-2+1] + wfnc[N-2-1] - 2.*wfnc[N-2])/(2.*A*A*m);

       //Gives error message if norm = inf or NaN
       if(res!=res || res==(comp) 0 )
              {
                     cout << "Problem with energy calculation";
              }
       tmp.clear();
       return res;
}

void record_snap(vector<comp>& wfnc, int num)//saves snapshots in /data/snapshots/
{
       char name [255];
       ofstream out;

       sprintf(name,"data/snapshots/wave_%d.txt",num);
       out.open(name, ios::out);
       out.precision(10);
       switch (SYMM)
       {
       case 1:
              for (int sx =0; sx<N;sx++)
       {
              out << sx << "       " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
       }
              break;
       case 2:
              for (int sx =0; sx<N;sx++)
       {
              out << -N/2. + sx << "     " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
       }
              break;
       default:
              break;
       }

       out.close();
}

void wfnc_out(vector <comp>& wfnc, char* label) //save wavefunction in /data
{
    ofstream out;
    sprintf(label,"data/wave_%s.txt",label);
    out.open(label, ios::out);
    out.precision(10);
    switch (SYMM)
    {
    case 1:
           for (int sx =0; sx<N;sx++)
    {
           out << sx << "       " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
    }
           break;
    case 2:
           for (int sx =0; sx<N;sx++)
    {
           out << -N/2. + sx << "     " << real(wfnc[sx]) << "   " << imag(wfnc[sx]) << "   " << real(conj(wfnc[sx])*wfnc[sx]) << endl;
    }
           break;
    default:
           break;
    }

    out.close();
}
//Algorithm for solving tridiagonal system
inline void LU_solve(comp coeff,vector<comp>& b, vector<comp>& d, vector<comp>& f)
{
  size_t n = d.size();
  off_diag_star[0] = -coeff / b[0];
  b_star[0] = d[0] / b[0];

  // Create the off_diag_star and b_star coefficients in the forward sweep
  for (size_t i=1; i<n; i++) {
    comp m = 1.0 / (b[i] - (-coeff) * off_diag_star[i-1]);
    off_diag_star[i] = -coeff * m;
    b_star[i] = (d[i] - (-coeff) * b_star[i-1]) * m;
  }

  // This is the reverse sweep, used to update the solution vector f
  for (int i=n-1; i-- > 0; ) {
    f[i] = b_star[i] - off_diag_star[i] * d[i+1];
  }
}

//set new b-vector
inline void update_b(comp a, comp t_imag,vector <comp>& wfnc)
{
       for (int i=1;i<N-1;i++)
       {
              b[i]=a*wfnc[i-1] + wfnc[i]*(1. - 0.5*t_imag*(1./(A*A*m)+V[i])) + a*wfnc[i+1];
       }
}

//Solving routine - Parameter a is for auto runs
void solve(double a)
{
       char name[255];comp lastenergy=1e9;double epsilon;int snaps=0;
       update_b(alpha,EPS,psi);
       ofstream out;
       ofstream dec ("data/decay.txt");
       if(RECORD==1){
       sprintf(name,"psi_t_%f",a);}
       else{
           sprintf(name,"ground_%f_%s",a,buf);
       }
       out.open(name, ios::out);
       out.precision(10);
       do
       {
           if (RECORD==1){ //Record decay if chosen to be
           if(t> 0.3*TIME && t<0.8*TIME)
           {
              dec << (double)(t*real(EPS)) << "\t" << real(psi[record_point]) <<"\t" << imag(psi[record_point]) << endl;
           }}
           else{//default: normalize
               energy=wfncEnergy(psi);

               if (t%snaptime==0){
                   normalize(psi);
                   energy=wfncEnergy(psi)/normalization;

                   epsilon =abs(abs(last-energy)/lastenergy);//check convergence
                   if (epsilon < 1e-8)
                   {
                         cout << "tolerance achieved!" << endl;
                         cout << "Difference: " << "\t" << epsilon << endl;
                         cout << "Energy: " <<"\t" << energy/normalization << endl;
                         break;
                   }
                   lastenergy=energy;
                   if (snaps < snap_num)
                   {
                      //record_snap(psi,snaps);
                      for(int i=0;i<N;i++)
                      {
                          psi_snap[i][snaps] = psi[i]; //save snap in data array
                      }
                      snaps++;

                   }

               }

           }
              LU_solve(alpha,diag,b,psi);
              update_b(alpha,EPS,psi);
              t++;
              cout << fixed;
              cout << setw(5) << "\r" << t << " " <<(double) t/TIME*100 <<" %  " << " " << "Psi("<<record_point<<"): " << real(psi[record_point]) <<flush;
       }
       while(t<=TIME);
       normalize(psi);
       energy0 = wfncEnergy(psi)/normalization;
       wfnc_out(psi,name);
       cout << endl;
       cout << "E_ground: " << energy0 << "\t alpha_s: " << alpha_s_im << endl;

       //save energies for different alpha_s in data file
       if(AUTO==1){
       ener1 << alpha_s_im << "\t" << real(energy0) << "\t" << imag(energy0) << endl;
       }
       update_time = clock() - update_time;
       out.close();
       dec.close();
}

void evolve(vector <comp>& init,vector<comp>& iter) //Evolve |PSI_1> = |newInit> - |newInit><psi_0|newInit>
{
    t=0;
    comp overlap;
    while(t<t_end){
    for (int i=0;i<N-1;i++)
      {
        overlap += A/2*(conj(init[i])*iter[i] + conj(init[i])*iter[i+1]);
      }
        overlap += -0.5*A*conj(init[N-1])*iter[N-1];

      for (int i=0;i<N;i++)
      {
          init[i] -= iter[i]*overlap;
      }
      if(t%snaptime==0){
               normalize(init);
      }
      LU_solve(alpha,diag,b,init);
      update_b(alpha,EPS,init);
      cout << fixed;
      cout << setw(5) << "\r" << t << "\t" <<(double) t/(t_end)*100 <<" %  "<<flush;
      t++;
      }
      normalize(init);
      energy1=wfncEnergy(init)/normalization;
}
void excited_state() //Routine to calculate 1st excited state
 {
     cout << endl;
     cout << "Calculating excited States" << endl;
     char name [255];
     ofstream out;
     for (int snaps=snap_num-1;snaps>snap_num-2;snaps--) //calculate 1st excited using last snapshot
     {
        cout << "Snapshot: " << snaps+1 << endl;
        sprintf(name,"data/excited/excited_%d.txt",snaps);
        out.open(name, ios::out);
        out.precision(10);
          for (int i=0;i<N;i++)
          {
              PSI_1[i]=psi_snap[i][snaps];
          }
        evolve(PSI_1,psi);

       for (int i=0;i<N;i++)
       {
             psi[i]=PSI_1[i];
       }
       normalize(PSI_1);
       energy1 = wfncEnergy(PSI_1)/normalization;

       for (int i=0;i<N-1;i++)
       {
             out << i*A << "\t" << real(PSI_1[i]) << "\t" << imag(PSI_1[i]) << "\t" << conj(PSI_1[i])*PSI_1[i] << endl;
       }
       out.close();
     }
}
void output() //Write info File
{
              cout << "Writing Info File" << endl;

              debug <<"Potential: " << buf << endl;
              debug <<"Initial: " << initbuf << endl;
              debug <<"Reduced Mass: " << m << endl;
              debug <<"Temperature: " << T << endl;
              debug << endl;
              debug << "-------------Results--------------- " << endl;
              debug << endl;
              debug <<"Energy of ground state: " << energy0 << endl;
              debug <<"Energy of 1st excited state: " << energy1 << endl;

              debug << endl;
              debug << "-------------Runtime--------------- " << endl;
              debug << endl;
              debug << "Iterations: " << t <<endl;
              debug << "Update Runtime: " << (float)update_time/(CLOCKS_PER_SEC*60)  << " min." << endl;
              debug << "Time per loop: " << (float)update_time/(CLOCKS_PER_SEC*t) << " sec."<< endl;
}

void automatize(){
    double p=1;//Vary values of alpha_s automatically
    while (p >0)
     {
         set_param();
         t=0;
         solve(p);
         output();
         if(p>0.2){
         p-=0.1;
         }
         else{
             p-=0.02;
         }
         alpha_s_im*=p;
         g = sqrt(4*M_PI*alpha_s_im);
         a_0 = 2/(m*C_f*alpha_s_im);
         m_D = pow(g*T,2)/3. *(N_c+(double)N_f/2.);
     }

}

int main()
{
       comp tmp;
       cout << setprecision(7);
       cout << "Choose potential: " << endl;
       cout << "(1) Coulomb " << endl;
       cout << "(2) Oscillator " << endl;
       cout << "(3) Particle in a box " << endl;
       cout << "(4) Quark - Antiquark " << endl;
       cout << ">> ";
       cin >> POTENTIAL;

       cout << "Choose initial wave: " << endl;
       cout << "(1) Gauss " << endl;
       cout << "(2) Delta peak " << endl;
       cout << "(3) Radial coulomb " <<endl;
       cout << "(4) Box ground state (sin()) " <<endl;
       cout << "(5) Advanced Dirac delta -6*N_c*r*d(r) " <<endl;
       cout << ">> ";
       cin >> INITCONDITION;

       cout << "Record psi(t) (1) or Renormalize (2-default)?" << endl;
       cout << ">>";
       cin >>RECORD;
       cout << "Automatize procedure (1)-(NO is default!)?" << endl;
       cout << ">>";
       cin >>AUTO;

       set_param();
       cout << endl;
       tmp = wfncEnergy(psi);

       //Output general Info and save them

       cout << " ------General Information------- " << endl;
       cout << endl;
       cout << "Gridpoints: " << N << " Spatial Step: " << A << " Time Step: " << EPS << " Alpha: " << real(alpha) << endl;
       cout << "Energy of initial wavefunction (without norm): " << tmp << endl;
       normalize(psi);
       cout << "Norm of inital wavefunction: " << normalization << endl;
       cout << "Normalized energy: " << tmp/normalization << endl;
       cout << endl;

       debug << " ------General Information------- " << endl;
       debug << endl;
       debug << "Gridpoints: " << N << " Spatial Step: " << A << " Time Step: " << EPS <<" Alpha: " << real(alpha) << endl;
       debug << "Record Point: " << record_point << endl;
       debug << "Energy of initial wavefunction (without norm): " << tmp << endl;
       debug<< "Norm of inital wavefunction: " << normalization << endl;
       debug << "Normalized energy: " << tmp/normalization << endl;
       debug << endl;

       cout << " ------Starting Crank Nicolson Evolution------- " << endl;
       cout << "Iterations" << "\t" << " % " << endl;

       if(AUTO==1){
       automatize();}
       else{
           set_param();
           solve(POTENTIAL);
           if (RECORD!=1){
           excited_state();}
           output();
       }

       cout << "****done***** " << endl;
       cin.get();
       return 0;
}
