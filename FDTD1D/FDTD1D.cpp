#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <cstring>
#include <cstdlib>
#include <stdlib.h>
using namespace std;

typedef complex<double> comp;
typedef vector<comp> cvec;
#ifndef M_PI
#define M_PI    3.1415926535897
#endif
//-----PHYSICAL CONSTANTS (QCD-SYS->c=1,hbar=1,4pi*e_0=1)---------//
const double m=4.14/2;                                        //MASS for evolution
const double hbar=1;                                          // h/2pi
const double e_0 =1;                                          //electric field constant
const double T = 0.4;                                        //Temperature
const short int N_c = 3;                                      //# colors
const short int N_f = 2;                                      //# quark flavours
const double C_f = (double)(N_c*N_c-1)/((double)(2*N_c));
const double alpha_s = 0.41;                                  //Strong coupling
double alpha_s_im = 0.41;                                     //strong coupling to manipulate
const double gamma_E = 0.57721566490153286060;                //euler constant
double g = sqrt(4*M_PI*alpha_s_im);
double a_0 = 2/(m*C_f*alpha_s);                               //'bohr' radius
double m_D = (pow(g*T,2)/3. *(N_c+(double)N_f/2.));           //Debye Mass

//------GRID and EVOLUTION-----------------------------
const size_t N=2000;                                          //Gridpoints
const double A = 0.05;                                        //Spatial Step
const double EPS = 0.0001;                                    //Time step
const int TIME =300000;                                      //Total time
int t;                                                        //Time
const int snap_num = 4;                                       //# of snapshots
const int snaptime = TIME/(snap_num);                        //time to check evolution and record snapshot
//------OTHER CONSTANTS AND VARS-----------------------------
comp im = comp(0,1);                                          //imaginary unit
char buf[255],initbuf[255];                                   //char to save info
comp normalization;                                           //holds norm
comp energy =1e9,energy1,energy0;                             //hold energies
int potoption,intoption,SYMM,RECORD,AUTO;                     //vars for options
clock_t update_time;                                          //check runtime
const int record_point =1;
//-----------FILES-------------------------
ofstream debug ("data/INFO.txt");
ofstream pot("data/potential.txt");
ofstream ener1("data/energies400.txt"); //files for qq-energies
//ofstream ener2("data/energies2.txt");
//----------Arrays and Matrices-------------
cvec a(N);                                                    //necessary fur update
cvec b(N);                                                    //-----------"---------
cvec psi(N);                                                  //holds the ground state
cvec V(N);                                                    //holds potential
vector < cvec > psi_snap(N,cvec(snap_num));                   //holds snapshots
cvec PSI_1(N);                                                //holds excited state


void wfnc_out(cvec& wfnc, char* label) //save wavefunction in /data
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

void set_pot()
{
    double r;
    //set pot parameters
    //set potential
    double angular = 0.;
    switch (potoption)
    {
    case 1: //coulomb
           SYMM = 1;
           sprintf(buf,"coulomb");
           for (int i=1;i<N;i++)
           {
               r=(double)i*A;
               V[i]=-pow(r,-1)+ angular*(angular+1)/(pow(r,2)*m*2.);
           }
           V[0]=(double)-1./(A*A);
           break;
    case 2://oscillator
           SYMM = 1;
           sprintf(buf,"osc");
           for (int i=0;i<N;i++)
           {
               r=((double)i-(double)N/2)*A;
               V[i]=0.5*pow(r,2);
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
    case 4://delta
           SYMM = 2;
           sprintf(buf,"delta");
           for (int i=0;i<N;i++)
           {
                  V[i]=0.;
           }
           V[N/2] = 1.;
           break;
    case 5: //QQ
              SYMM = 1;
              sprintf(buf,"QQ");
              for (int i=1;i<N;i++)
              {
                  r=i*A;
                  V[i]=(comp)-C_f*alpha_s*pow(r,-1) +(comp)im*((C_f*(alpha_s_im*pow(r,2)*T*pow(m_D,2)))/((comp)6.) * (2*gamma_E+ log(pow(r*m_D,2))-(comp)8./3.));
              }
              V[0]=-1./(pow(A,2));
              break;
    default:
           break;
    }
    for (int i=0;i<N;i++)
    {
        pot << i*A << "\t" << real(V[i]) << "\t" << imag(V[i]) << endl;
    }
    for(int i=0;i<N;i++)
    {
        a[i] = (1.-EPS/2.*V[i])/(1.+EPS/2.*V[i]);
        b[i] = 1./(1.+EPS/2.*V[i]);

    }



}

void set_init()
{
    double x,sgn;
    //-------set initial wave---------------

switch (intoption)
{
case 1: //random

    sprintf(initbuf,"Random");
    for (int i=0;i<N;i++)
    {
        sgn=(double)1./(rand()%100 +1);
        if(sgn>0.5)
        {
           psi[i]=(double)-1./(rand()%100 +1);                  }
        else{
            psi[i]=(double)1./(rand()%100 +1);
        }
    }
    break;
case 2: //delta peak
    sprintf(initbuf,"delta peak");
    for (int i=0;i<N;i++)
    {
           if(i==N/2)
           {
                 psi[N/2]=1;
           }
           else
           {
                 psi[i]=0.0;
           }
    }
    break;
case 3://radial of coulomb
    sprintf(initbuf,"radial coulomb");
    for (int i=0;i<N;i++)
    {
           x=i*A;
           psi[i]=2*x*exp(-x);
    }
    break;
case 4: //Particle in a box ground state
    sprintf(initbuf,"box ground state");
    for (int i=0;i<N;i++)
    {
           x=i*A;
           psi[i]= sqrt(2./(A*N)) * sin(2*M_PI * x /(A*N));
    }
    break;
case 5: //advanced dirac delta
    sprintf(initbuf,"advanced dirac delta");
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

void set_param()//Set conditions
{
set_pot();
set_init();
}

inline void normalize(cvec& wfnc)//normalizes wavefunction
{
       comp sum;
       //via trapez rule
       cvec tmp (2);
       for (int i=0;i<N-1;i++)
       {
              //cout << i << " " << wfnc[i] << endl;
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
       if (sum!=sum || sum==(comp)0)
       {
              cout << "\a Problem with normalization!";
       }
       tmp.clear();
}

comp wfncEnergy(cvec& wfnc)//computes energy
{
       complex <double> res=0;
       //via trapez
       cvec tmp(2);
       for (int sx=1;sx<N-2;sx++)
              {
                     tmp[0] = V[sx]*conj(wfnc[sx])*wfnc[sx]- conj(wfnc[sx])*hbar*hbar*(wfnc[sx+1] + wfnc[sx-1] - 2.*wfnc[sx])/(2.*A*A*m);
                     tmp[1] = V[sx+1]*conj(wfnc[sx+1])*wfnc[sx+1]- conj(wfnc[sx+1])*hbar*hbar*(wfnc[sx+1+1] + wfnc[sx-1+1] - 2.*wfnc[sx+1])/(2.*A*A*m);
                  res += A*(tmp[0]+0.5*(tmp[1]-tmp[0]));
              }
       res+= -0.5*A*V[N-2]*conj(wfnc[N-2])*wfnc[N-2]- conj(wfnc[N-2])*hbar*hbar*(wfnc[N-2+1] + wfnc[N-2-1] - 2.*wfnc[N-2])/(2.*A*A*m);
       if(res!=res || res==(comp) 0 )
              {
                     cout << "\a Problem with energy calculation";
              }
       tmp.clear();
       return res;
}

void record_snap(cvec& wfnc, int num)//saves snapshots in /data/snapshots/
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


inline void update(cvec &wfnc) //Update routine
{
    for (int i=1;i<N;i++)
    {
        wfnc[i] = a[i]*wfnc[i] + b[i]*EPS/(2.*A*A*m)*(wfnc[i-1]+wfnc[i+1]-2.*wfnc[i]);
    }
}

inline void update_boundary(cvec &wfnc)//Incorporate Boundaries
{
 wfnc[0]=0;
 wfnc[N-1]=0;
}

void solve(double a) //Solving routine
{
       int snaps = 0;comp temp_energy = 1e9;double epsilon=1e9;char name[255];
       update_time=clock();
       ofstream out;
       ofstream dec ("data/decay.txt");
       if(RECORD==1){
       sprintf(name,"data/psi_t_%f.txt",a);}
       else{
           sprintf(name,"ground_%f_%s.txt",a,buf);
       }

       t=1;
       do
       {
           if (RECORD==1){ //Record decay if chosen to be
           if(t> 0.3*TIME && t<0.8*TIME)
           {
              dec << (double)(t*EPS) << "\t" << real(psi[record_point]) <<"\t" << imag(psi[record_point]) << endl;
           }}
           else{//default:normalize
              energy = wfncEnergy(psi);
              if (t%snaptime==0)
              {
                     normalize(psi);
                     epsilon =abs(abs(temp_energy-energy)/temp_energy);
                     if (epsilon < 1e-8)
                     {
                           cout << "tolerance achieved!" << endl;
                           cout << "Difference: " << "\t" << epsilon << endl;
                           cout << "Energy: " <<"\t" << energy/normalization << endl;
                           break;
                     }

                     temp_energy=energy;
                     if (snaps < snap_num)
                     {
                         //record_snap(psi,snaps);
                        for(int i=0;i<N;i++)
                        {
                            psi_snap[i][snaps] = psi[i];
                        }
                        snaps++;
                     }
               }}
              update(psi);
              update_boundary(psi);
              cout << fixed;
              cout << setw(5) << "\r" << t << "\t   \t" <<(double) t/TIME*100 <<" %  "<<flush;
              t++;
       }
       while(t<TIME);
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

void evolve(cvec& init,cvec& iter) //Evolve |PSI_1> = |newInit> - |newInit><psi_0|newInit>
{
    t=0;
    comp overlap;
    while(t<TIME){
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
      update(init);
      update_boundary(init);
      cout << fixed;
      cout << setw(5) << "\r" << t << "\t" <<(double) t/(TIME)*100 <<" %  "<<flush;
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
     //ener2 << alpha_s_im <<"\t" << real(energy1)<<"\t" <<imag(energy1) <<endl;


 }

void output() //Write information File
{
              cout << "Writing Info File" << endl;

              debug <<"Potential: " << buf << endl;
              debug <<"Initial: " << initbuf << endl;
              debug <<"Mass: " << m << endl;
              debug <<"Temperature: " << T << endl;
              debug << endl;
              debug << "-------------Results--------------- " << endl;
              debug << endl;
              debug <<"Energy of ground state: " << energy << endl;
              debug <<"Energy of excited state: "<< energy1 << endl;

              debug << endl;
              debug << "-------------Runtime--------------- " << endl;
              debug << endl;
              debug << "Iterations: " << t <<endl;
              debug << "Update Runtime: " << (float)update_time/(CLOCKS_PER_SEC*60)  << " min." << endl;
              debug << "Time per loop: " << (float)update_time/(CLOCKS_PER_SEC*t) << " sec."<< endl;
}

void automatize(){

    double p=1.;  //automatize procedure
    while(p>0)
    {
        if(p>0.1){
            p-=0.05;
        }
        else{
            p-=0.01;
        }
        alpha_s_im =0.41*p;
         g = sqrt(4*M_PI*alpha_s_im);
         //a_0 = 2/(m*C_f*alpha_s_im);
         m_D = pow(g*T,2)/3. *(N_c+(double)N_f/2.);
         set_param();
         solve(p);
         //excited_state();
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
       cout << "(4) Delta" << endl;
       cout << "(5) Quark-Antiquark " << endl;
       cout << ">> ";
       cin >> potoption;

       cout << "Choose initial wave: " << endl;
       cout << "(1) Random " << endl;
       cout << "(2) Delta peak " << endl;
       cout << "(3) Radial coulomb " <<endl;
       cout << "(4) Box ground state (peak sin()) " <<endl;
       cout << "(5) Advanced Dirac delta -6*N_c*r*d(r) " <<endl;
       cout << ">> ";
       cin >> intoption;
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
       cout << "Gridpoints: " << N << " Spatial Step: " << A << " Time Step: " << EPS << endl;
       cout << "Energy of initial wavefunction (without norm): " << tmp << endl;
       normalize(psi);
       cout << "Norm of inital wavefunction: " << normalization << endl;
       cout << "Normalized energy: " << tmp/normalization << endl;
       cout << endl;

       debug << " ------General Information------- " << endl;
       debug << endl;
       debug << "Gridpoints: " << N << " Spatial Step: " << A << " Time Step: " << EPS << endl;
       debug << "Energy of initial wavefunction (without norm): " << tmp << endl;
       debug<< "Norm of inital wavefunction: " << normalization << endl;
       debug << "Normalized energy: " << tmp/normalization << endl;
       debug << endl;

       cout << " ------Starting Explicit Evolution------- " << endl;
       cout << "Iterations" << "\t" << " % " << endl;

       if(AUTO==1){
       automatize();}
       else{
           set_param();
           solve(potoption);
           if(RECORD=!1){
           excited_state();}
           output();
       }

       cout << "****done***** " << endl;
       cin.get();
       return 0;
}
