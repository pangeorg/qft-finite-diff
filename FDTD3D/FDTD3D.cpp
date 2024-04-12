#include <iostream>
#include <cstdlib>
#include <time.h>
#include <complex>
#include <fstream>
#include <iomanip>

using namespace std;

typedef complex<double> comp;
#ifndef M_PI
#define M_PI    3.1415926535897
#endif
//--------PHYSICAL CONSTANTS--------------------
const double omega = 1.;
const double hbar=1;//1.05457148e-34;                         // h/2M_PI
const double ec= 1;//-1.60221e-19;                            //charge
double m=1;                                                   //Mass used for calculation
const double e_0 =1;//8.854187817e-12;                        //electric field constant
const double T = 782.778;                                           //Temperature
const double beta=1/T;                                        // 1/T
const short int N_c = 3;
const short int N_f = 2;                                      //# quark flavours
const double C_f = 1;
const double alpha_s = 0.13;                                     //Strong coupling
const double gamma_E = 1;
const double g = sqrt(4*M_PI*alpha_s);
const double a_0 = 1;
const double m_e = 1;//9.10938188e-31;                        //Electron Mass
const double m_D = pow(g*T,2)/3. *(N_c+(double)N_f/2.);
const comp im = (comp)(0,1);
//-------------GRID-----------------------
const unsigned int NUM=170;				//Number of Spatial Steps
const unsigned int TIME= 40000;			//Time Steps
const unsigned int snaptime=1000;		//Snapshot #
const double A=0.05;					//Spatial step
const double EPS=0.001;//A*A/8;         //time step should be smaller than A^2/3 to achieve stability
unsigned int t = 0;						//imaginary time
//-----------USEFUL----------------------
unsigned int counter=0;
int potoption,initoption;
const int UPDATE = 100;
//nsnaps = t/SNAPTIME;
int nsnaps = 0;
int snaps=nsnaps+1;
const int INITSYMMETRY=1;
double FDTDTime;
double tauf = T*EPS;					//tau >>1/(E_2-E_1)
double TOLERANCE = 1e-9;
double convtoeV =1;// 27.2;				//conversion from at to eV
double Vmax = 1.;
clock_t progtimeFDTD;
comp energytot,lastenergy = 1e10;
comp normalization,waveenergy,over,energy_1;
//-------------WAVES----------------------
comp*** psi;                            //current wavefunktion
comp*** PSI;                            //updated wavefunktion
comp*** PSI_1;                          //1st exited state
comp*** a;                              //coefficents in update rule
comp*** b;
comp**** psisnap;                       //snapshots of wavefunction
comp*** v;                              //Potential array FEPSD
comp*** tmp;
comp* energies;                         //energies of higher states

//-----------FILES-------------------//
ofstream debug ("debug.txt");

//--------Functions-----------------//



double gauss(double x,double t,double sigma)
{

    return exp(-((x-t)*(x-t))/(2*sigma*sigma));
}
//--------POTENTIAL-----------------//
comp V(double x,double y,double z)
{
    double r;
    double xx,yy,zz;

    xx = ((double) x) - ((double)NUM)/2.;
    yy = ((double) y) - ((double)NUM)/2.;
    zz = ((double) z) - ((double)NUM)/2.;
    r = A*sqrt(xx*xx+yy*yy+zz*zz);
    switch (potoption)
    {
    case 1:
        return -Vmax;								//Finite Potential well
        break;
    case 2:
            return 0.5*m*omega*omega*xx*xx*A*A;	//Oscillator
        break;
    case 3:
        return 0.;									//Infinit Potential well
        break;
    case 4:
        if (r == 0)
            return -1./(pow(A,2));
        else
            return (comp)-pow(r,-1) + (comp)im*((C_f*(alpha_s*pow(r,2)*T*pow(m_D,2)))/((comp)6.) * (2*gamma_E+pow(log(r*m_D),2)-(comp)8./3.));						//Quarkonium Potential
        break;
    case 5:											
		if(r==0)
		{
            return -1./(pow(A,2));
		}
		else
		{
        return -(comp)1./r;								//coulomb
		}
        break;
    default:
        break;
    }
}

void load_potential()
{
        int sx,sy,sz;
        for (sx=0;sx<=NUM;sx++)
        for (sy=0;sy<=NUM;sy++)
        for (sz=0; sz<=NUM;sz++) {
            v[sx][sy][sz] = V(sx,sy,sz);
                b[sx][sy][sz] = 1./(1.+EPS*v[sx][sy][sz]/((comp) 2.));
                a[sx][sy][sz] = (1.-EPS*v[sx][sy][sz]/((comp) 2.))*b[sx][sy][sz];
        }
}
void init_wfkt()
{
    double r,x,y,z;

    switch (initoption)
    {
    case 1:
        cout << "Constant " << endl;
        for (int i=0;i<=NUM;i++)
            for (int j=0;j<=NUM;j++)
                for (int k=0;k<=NUM;k++)
                {
                    psi[i][j][k] = 0.2;
                }
        break;
    case 2:
        cout << "gaussian wave " << endl;
        for (int i=0;i<NUM+1;i++)
            for (int j=0;j<NUM+1;j++)
                for (int k=0;k<NUM+1;k++)
                {
                    r = A*sqrt(i*i+j*j+k*k);
                    psi[i][j][k] = (comp) gauss(r,NUM/2,4);
                }
        break;
    case 3:
        cout << "Sinus " << endl;
        for (int i=0;i<NUM+1;i++)
            for (int j=0;j<NUM+1;j++)
                for (int k=0;k<NUM+1;k++)
                {
                    psi[i][j][k] = sin((double)(M_PI*i*1/NUM));
                }
        break;
    case 4:
        cout << "delta peak " << endl;
        for (int i=0;i<NUM+1;i++)
            for (int j=0;j<NUM+1;j++)
                for (int k=0;k<NUM+1;k++)
                {
                    psi[i][j][k] = 0;
                }
        psi[NUM/2][NUM/2][NUM/2]=1;
        break;
    case 5:
        cout << "Coulomb, l,m=0" << endl;
        for (int i=0;i<NUM+1;i++)
            for (int j=0;j<NUM+1;j++)
                for (int k=0;k<NUM+1;k++)
                {
                    x = (double)i - (double)(NUM/2.);
                    y = (double)j - (double)(NUM/2.);
                    z = (double)k - (double)(NUM/2.);
                    r = A*sqrt(x*x+y*y+z*z);
                    psi[i][j][k] = 2.*exp(-r)*(comp)1./(sqrt(4*M_PI));
                }
        break;
        default:
        break;
    }
    int sx,sy,sz;
    // enforce BCs
        for (sx=0;sx<=NUM;sx++)
          for (sy=0;sy<=NUM;sy++) {
        psi[sx][sy][0] = 0;
        psi[sx][sy][NUM] = 0;
      }

        for (sz=0;sz<=NUM;sz++)
          for (sy=0;sy<=NUM;sy++) {
        psi[0][sy][sz] = 0;
        psi[NUM][sy][sz] = 0;
      }

        for (sx=0;sx<=NUM;sx++)
          for (sz=0;sz<=NUM;sz++) {
        psi[sx][0][sz] = 0;
        psi[sx][NUM][sz] = 0;
      }


    // zero out updated wavefnc for safety's sake
        for (int sx=0;sx<NUM;sx++)
          for (int sy=0;sy<NUM;sy++)
            for (int sz=0;sz<NUM;sz++)
                PSI[sx][sy][sz] = 0;

}
void allocate_memory()
{
    cout << "Preparing Memory..." << endl;

    psi = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) psi[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) psi[sx][sy] = new comp[NUM];

    PSI = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) PSI[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) PSI[sx][sy] = new comp[NUM];

    PSI_1 = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) PSI_1[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) PSI_1[sx][sy] = new comp[NUM];

    v = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) v[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) v[sx][sy] = new comp[NUM];

    a = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) a[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) a[sx][sy] = new comp[NUM];

    b = new comp**[NUM];
    for (int sx=0;sx<NUM+1;sx++) b[sx] = new comp*[NUM];
    for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) b[sx][sy] = new comp[NUM];

    //nsnaps = STEPS/SNAPUPDATE;
    int nsnaps = 0;
    int snaps = nsnaps+1;
    psisnap = new comp***[snaps];
    for (int n=0;n<=snaps;n++) psisnap[n] = new comp**[NUM];
    for (int n=0;n<=snaps;n++) for (int sx=0;sx<NUM+1;sx++) psisnap[n][sx] = new comp*[NUM];
    for (int n=0;n<=snaps;n++) for (int sx=0;sx<NUM+1;sx++) for (int sy=0;sy<NUM+1;sy++) psisnap[n][sx][sy] = new comp[NUM];

    energies = new comp [snaps];

    cout << "Done" << endl;
}

//normalize wavefunction
void normalize(comp*** f,int option) {
    comp norm;
    switch(option)
    {
    case 1://trapez
        for (int sx=0;sx<NUM;sx++)
          for (int sy=0;sy<NUM;sy++)
            for (int sz=0; sz<NUM;sz++) {
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x,y,z)
                sx +=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x+1,y,z)
                sx-=1;sy+=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x,y+1,z)
                sy-=1;sz+=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x,y,z+1)
                sz-=1;sx+=1;sy+=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x+1,y+1,z)
                sz+=1;sy-=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x+1,y,z+1)
                sx-=1;sy+=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x,y+1,z+1)
                sx+=1;
                norm+=(conj(f[sx][sy][sz]))*f[sx][sy][sz]; //f(x+1,y+1,z+1)
                sx-=1;sy-=1;sz-=1;

            }
            norm*=pow(A,3)/((double)8.);
            break;
    case 2:
        for (int sx=0;sx<=NUM;sx++)
          for (int sy=0;sy<=NUM;sy++)
            for (int sz=0; sz<=NUM;sz++) {
                norm +=f[sx][sy][sz];
            }
            norm*=pow(A,3)/((double)8.);
        break;
     default:
        break;

    }

    comp norm_sqrt = sqrt(norm);
        for (int sx=0;sx<=NUM;sx++)
          for (int sy=0;sy<=NUM;sy++)
            for (int sz=0; sz<=NUM;sz++) {
                f[sx][sy][sz]/=norm_sqrt;
            }
        normalization=norm;
        if (norm == (comp) 0)
        {
            cout << "Norm = 0 ! " << endl;
        }
        if (norm!=norm)
        {
            cout << "QNAN norm" << endl;
        }
}
void recordsnap(comp*** wfnc, int step) {
    int snap = (int) step/snaptime;
        for (int sx=0;sx<=NUM;sx++)
          for (int sy=0;sy<=NUM;sy++)
            for (int sz=0; sz<=NUM;sz++) {
                psisnap[snap][sx][sy][sz] = wfnc[sx][sy][sz];}
}
inline comp update(int sx, int sy, int sz, double step)
{
  return psi[sx][sy][sz]*a[sx][sy][sz] + b[sx][sy][sz]*step*hbar*hbar*(
          psi[sx+1][sy][sz] + psi[sx-1][sy][sz] +
          psi[sx][sy+1][sz] + psi[sx][sy-1][sz] +
          psi[sx][sy][sz+1] + psi[sx][sy][sz-1] -
          ((comp) 6.)*psi[sx][sy][sz])/(((comp) 2.)*A*A*m);
}

//-----compute energy of wavefunction-------//
inline comp wfncEnergy(comp*** wfnc, int option) {
    comp res=0;
    double coeff = pow(hbar,2)/(2.*m*A*A);
    switch(option)
    {
    case 1:    //via trapez
        for (int sx=1;sx<NUM-2;sx++)
          for (int sy=1;sy<NUM-2;sy++)
            for (int sz=1;sz<NUM-2;sz++)
            {
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x,y,z)
                sx +=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x+1,y,z)
                sx-=1;sy+=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x,y+1,z)
                sy-=1;sz+=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x,y,z+1)
                sz-=1;sx+=1;sy+=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x+1,y+1,z)
                sz+=1;sy-=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x+1,y,z+1)
                sx-=1;sy+=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x,y+1,z+1)
                sx+=1;
                res+=(conj(wfnc[sx][sy][sz]))*(v[sx][sy][sz] * wfnc[sx][sy][sz] - coeff  * (wfnc[sx+1][sy][sz]+wfnc[sx-1][sy][sz]+wfnc[sx][sy+1][sz]+wfnc[sx][sy-1][sz]+wfnc[sx][sy][sz+1]+wfnc[sx][sy][sz-1] - 6.*wfnc[sx][sy][sz])); //f(x+1,y+1,z+1)
                sx-=1;sy-=1;sz-=1;

            }
            res*=pow(A,3)/8.;
        break;
    case 2: //sum (fast)

        for (int sx=1;sx<NUM;sx++)
          for (int sy=1;sy<NUM;sy++)
            for (int sz=1;sz<NUM;sz++)
            {
                res += v[sx][sy][sz]*conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz] -
                   conj(wfnc[sx][sy][sz])*(  wfnc[sx+1][sy][sz] + wfnc[sx-1][sy][sz] +
                                     wfnc[sx][sy+1][sz] + wfnc[sx][sy-1][sz] +
                                     wfnc[sx][sy][sz+1] + wfnc[sx][sy][sz-1] -
                                     ((comp) 6.)*wfnc[sx][sy][sz] )*coeff;
            }
        res*=pow(A,3)/((double)8.);
        break;
    default:
        break;
    }
    return res;
}

inline void copyDown() {
    tmp = psi;
    psi = PSI;
    PSI = tmp;
}
void outputSnapshot(comp ***wfnc, char* label) {

  double x;
  static int h=NUM/2;
  static int hx=NUM/2;

  ofstream out;
  char fname[255];

  // dump wavefunction

  // output radial slice suitable for 2d viewing
  switch(potoption)
  {
      case 1:
          sprintf(fname,"data/snapshot/wavefunction_pot_well_%s.txt",label);
          break;
      case 2:
          sprintf(fname,"data/snapshot/wavefunction_osc_%s.txt",label);
          break;
      case 3:
          sprintf(fname,"data/snapshot/wavefunction_free_%s.txt",label);
          break;
      case 4:
          sprintf(fname,"data/snapshot/wavefunction_cornell_%s.txt",label);
          break;
      case 5:
          sprintf(fname,"data/snapshot/wavefunction_coulomb_%s.txt",label);
          break;
      default:
          break;
    }
  out.open(fname, ios::out);
  out.precision(10);
  for (int s=0;s<=NUM;s++) {
    x=A*s;
    out << x<< "\t";
    out << scientific << real (0.5*(wfnc[s][h][h]+wfnc[s][h+1][h+1])) << " ";
    out << scientific << imag (0.5*(wfnc[s][h][h]+wfnc[s][h+1][h+1])) << "\t";
    out << endl;
  }
  out << "&&" << endl;
  for (int s=0;s<=NUM;s++) {
      x=A*s;
      out << x<< "\t";
    out << scientific <<real( 0.5*(wfnc[hx][s][h]+wfnc[hx+1][s][h+1])) << " ";
    out << scientific <<imag( 0.5*(wfnc[hx][s][h]+wfnc[hx+1][s][h+1])) << "\t";
    out << endl;
  }    over*=pow(A,3)/8.;
  out << "&&" << endl;
  for (int s=0;s<=NUM;s++) {
      x=A*s;
      out << x<< "\t";
    out << scientific <<real( 0.5*(wfnc[hx][h][s]+wfnc[hx+1][h+1][s])) << " ";
    out << scientific <<imag( 0.5*(wfnc[hx][h][s]+wfnc[hx+1][h+1][s])) << "\t";
    out << endl;
  }
  out.close();

}
void outputWavefunction(comp ***wfnc, char* label) {

  double x,y,z;
  ofstream out;
  char fname[255];
  // output full 3d wfnc
    switch(potoption)
  {
      case 1:
          sprintf(fname,"data/wavefunction_pot_well_%s.txt",label);
          break;
      case 2:
          sprintf(fname,"data/wavefunction_osc_%s.txt",label);
          break;
      case 3:
          sprintf(fname,"data/wavefunction_free_%s.txt",label);
          break;
      case 4:
          sprintf(fname,"data/wavefunction_cornell_%s.txt",label);
          break;
      case 5:
          sprintf(fname,"data/wavefunction_coulomb_%s.txt",label);
          break;
      default:
          break;
    }

  cout << "==> Dumping wave function to " << fname << endl;

  out.open(fname, ios::out);
  out.precision(12);
  for (int sx=0;sx<=NUM;sx++) {
    x=A*sx;
    for (int sy=0;sy<=NUM;sy++) {
        y=A*sy;
      for (int sz=0; sz<=NUM;sz++) {
          z=A*sz;
                out << x << "\t";
                out << y << "\t";
                out << z << "\t";
                out << real(wfnc[sx][sy][sz]) << "\t";
                out << imag(wfnc[sx][sy][sz]) << "\t";
                out << real(conj(wfnc[sx][sy][sz])*wfnc[sx][sy][sz]) << "\t";
                out << endl;
  }}}
  cout << "done" << endl;
  out.close();
}

void findExcitedStates() {
    char label[255];
    sprintf(label,"0_");
    outputWavefunction(psi,label);

    cout << "Calculating Excited state...." << endl;
    // compute first excited state
    for (int snap = 0;snap < snaps;snap++){ //int snap = nsnaps-1;
    // compute overlap
    //via trapez
    for (int sx=1;sx<NUM-1;sx++)
      for (int sy=1;sy<NUM-1;sy++)
        for (int sz=1;sz<NUM-1;sz++)
        {
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x,y,z)
            sx +=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x+1,y,z)
            sx-=1;sy+=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x,y+1,z)
            sy-=1;sz+=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x,y,z+1)
            sz-=1;sx+=1;sy+=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x+1,y+1,z)
            sz+=1;sy-=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x+1,y,z+1)
            sx-=1;sy+=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x,y+1,z+1)
            sx+=1;
            over+= psi[sx][sy][sz]*psisnap[snap][sx][sy][sz]; //f(x+1,y+1,z+1)
            sx-=1;sy-=1;sz-=1;

        }
    over*=pow(A,3)/((double)8.);

    // subtract overlap
    for (int sx=0;sx<=NUM;sx++)
        for (int sy=0;sy<=NUM;sy++)
            for (int sz=0; sz<=NUM;sz++){
                PSI_1[sx][sy][sz] = psisnap[snap][sx][sy][sz] - over*psi[sx][sy][sz];}

    normalize(PSI_1,1);
    energies[snap]=wfncEnergy(PSI_1,1)/normalization;

    sprintf(label,"1_%f",real(energies[snap]));
    outputWavefunction(PSI_1,label);
    }
}
void FDTD()
{
    progtimeFDTD =clock();
    char label [64];
    sprintf(label,"initial");
    cout << fixed;
    cout << endl;
    cout << "-------Starting Time Evolution---------" << endl;
    cout << endl;
    outputWavefunction(psi,label);
    do
    {
        waveenergy=wfncEnergy(psi,1)/normalization;
        if (t%snaptime==0)
        {
            normalize(psi,1);

            sprintf(label,"%d",t);
            //outputSnapshot(psi,label);
            energytot=(comp)waveenergy/normalization;

                cout << "Energy: " << energytot << endl;

            if (abs(energytot-lastenergy)<TOLERANCE)
            {
                cout << endl;
                cout << "t%snaptime =0 " << endl;
                cout << "Time: " << EPS*t << endl;
                cout << "Energy: " <<"\f" << energytot<< endl;
                cout << "\a";
                break;
            }
            else
            {
                lastenergy=energytot;
                if (t!=T)
                {
                    recordsnap(psi,0);
                }
            }
        }


        if (t<=TIME)
        {
            for(int i = 0; i<=UPDATE;i++){

               for (int sx=1;sx<NUM;sx++){
                  for (int sy=1;sy<NUM;sy++){
                    for (int sz=1;sz<NUM;sz++) {
                        PSI[sx][sy][sz] = update(sx,sy,sz,EPS);
                    }
                  }
                }
			}
             tmp=psi; //copy old wavefuntion
             psi=PSI;
             PSI=tmp;
         
        }
        cout << "\r"<< (double)t/TIME *100 << " % \t"<<flush ;
        counter++;
        t+=UPDATE;
    }
    while (t <= TIME);

    progtimeFDTD=clock()-progtimeFDTD;
    FDTDTime = (float) progtimeFDTD/CLOCKS_PER_SEC *1/60;
    normalize(psi,1);
    energy_1 = wfncEnergy(psi,1)/normalization;
    findExcitedStates();
}
void output()
{

    cout << "writing data..." << endl;

    debug <<"-------------General Information-------------" << endl;
    debug << endl;
    debug << "Gridpoints: " << NUM <<"^3 " << endl;
    debug << "Spatial Step: " << A << endl;
    debug << "Time Step: " << EPS << endl;

    switch (initoption)
    {
    case 1:
        debug << "Initialized wavefunction with constant" << endl;
        break;
    case 2:
        debug << "Initialized wavefunction with Gauss-Function" << endl;
        break;
    case 3:
        debug << "Initialized wavefunction with Sinus" << endl;
        break;
    case 4:
        debug << "Initialized wavefunction with Delta Peak" << endl;
        break;
    case 5:
        debug << "Initialized wavefunction with Coulomb" << endl;
        break;
    default:
        break;
    }

    switch (potoption)
    {
    case 1:
        debug << "Potential: Finite Potential well" << endl;
        break;
    case 2:
        debug << "Potential: Oscillator" << endl;
        break;
    case 3:
        debug << "Potential: Infinit Potential well" << endl;
        break;
    case 4:
        debug << "Potential: Quarkonium" << endl;
        break;
    case 5:
        debug << "Potential: Coulomb" << endl;
        break;
    default:
        break;
    }

    debug << "---------FDTD---------" << endl;

    debug << "OVERLAP: " << over << endl;

    debug <<"Energy E_0: " << real(energy_1) << " + i* " << imag(energy_1)<< endl;

    for (int i=0;i<=snaps;i++)
    {
        debug << "Energy_" << i+1 << ": " << energies[i] << endl;
        debug << "(E_" << i+1 << " - E_0)*tauf = " << (energies[i]-energy_1)*tauf << endl;

        if (real((energies[i]-energy_1)*tauf ) <4.)
        {
            debug << "tau_f too small, states degenerate" << endl;
        }
    }


    debug << "Runtime while Updating: " << FDTDTime << " min." << endl;
    debug << "Time per loop: " << FDTDTime *1/counter << " sec." << endl;


}
void input()
{
    cout << "----------Schrodinger Equation Solver with FDTD----------" << endl;
    cout << "Choose Potential: " << endl;
    cout << "(1) Finite Potential Well " << endl;
    cout << "(2) Harmonic Oscillator " << endl;
    cout << "(3) Free Particle " << endl;
    cout << "(4) Quarkonium Potential " << endl;
    cout << "(5) Coulomb Potential " << endl;
    cin >> potoption;

    cout << "Choose Initial Wave: " << endl;
    cout << "(1) Constant of 0.2 " << endl;
    cout << "(2) Gaussian Curve " << endl;
    cout << "(3) Sinus-Box " << endl;
    cout << "(4) Delta Peak " << endl;
    cout << "(5) Coulomb wave " << endl;
    cin >> initoption;
    cout << endl;
}

int main()
{

    cout << setprecision(7);
    comp tmp;
    input();
    allocate_memory();
    load_potential();
    init_wfkt();
    tmp=wfncEnergy(psi,1);
    cout << "--------Starting Information--------" << endl;
    cout << "Gridpoints in 3 Dimensions: " << NUM << " Spatial Step: " << A << " Time Step: " << EPS << endl;
    cout << "Energy of inital wavefunction: " << tmp << endl;
    normalize(psi,1);
    cout << "Norm of inital wave: " << normalization << endl;
    cout << "Normalized Energy: " << tmp/normalization << endl;
    FDTD();

    output();
    /*
    potential.close();
    wave.close();
    initwave.close();
    snapshot.close();
    */
    debug.close();

    cout << "----FINISHED------"<< endl;
    cin.get();
    return 0;
}
