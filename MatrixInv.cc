#include <iostream>
#include <cmath>
#define EPS 0.000000001
#define MAXIT 200000
using namespace std;
/* Using Jacobi iteration method */
void Jacobi(int n, float a[4][4], float b[4], float x[4], int *count, int *status){
  int    key;
  float  sum, x0[4];
  for(int i=0;i<4;i++)
    x0[i] = b[i] / a[i][i];

  *count = 1;
  begin:
  key = 0;

  for(int i=0;i<4;i++){
    sum = b[i];
    for(int j=0;j<4;j++){
      if(i==j)continue;
      sum = sum - a[i][j] * x0[j];
    }
    x[i] = sum / a[i][i];
    if(key == 0){
      /* Testing for accuracy */
      if(std::abs ((x[i]-x0[i]) / x[i]) > EPS ) key = 1;
    }
  }
  /* Testing for convergence */ 
  if(key == 1){
    if(*count == MAXIT){
      *status = 2;
      return; 
    }
    else{
      *status = 1;
      for (int i=0;i<4;i++)
        x0[i] = x[i];
    }
    *count = *count+1;
    goto begin;
  }
  return;

}
int main(){
  int     Npp, Npf,Nfp, Nff;
  float   ep1, ep2, f1, f2;
  float   Nyy, Nyj, Njy, Njj;
  float   a[4][4], b[4], x[4];
  int     n, count, status;  
  

  /* getting the input values */
  cout << "How many linear equations?"; cin >> n;
  cout << "Enter the values:" << endl;
  cout << "Npp: "; cin >> Npp;
  cout << "Npf: "; cin >> Npf;
  cout << "Nfp: "; cin >> Nfp;
  cout << "Nff: "; cin >> Nff;
  cout << "ep1: "; cin >> ep1;
  cout << "ep2: "; cin >> ep2;
  cout << "f1: ";  cin >> f1;
  cout << "f2: ";  cin >> f2;

  cout << "The equation to solve:" << endl;
  cout << Npp << "=" << ep1*ep2 << "*Nyy+" << ep1*f2 << "*Nyj+" << f1*ep2 << "*Njy+" << f1*f2 << "*Njj" << endl;
  cout << Npf << "=" << ep1*(1-ep2) << "*Nyy+" << ep1*(1-f2) << "*Nyj+" << (1-f1)*ep2 << "*Njy+" << f1*(1-f2) << "*Njj" << endl;
  cout << Nfp << "=" << (1-ep1)*ep2 << "*Nyy+" << (1-ep1)*f2 << "*Nyj+" << (1-f1)*ep2 << "*Njy+" << (1-f1)*f2 << "*Njj" << endl;
  cout << Nff << "=" << (1-ep1)*(1-ep2) << "*Nyy+" << (1-ep1)*(1-f2) << "*Nyj+" << (1-f1)*(1-ep2) << "*Njy+" << (1-f1)*(1-f2) << "*Njj" << endl;
  

  /* putting the values inside the arrays , remember ax=b */
  a[0][0] = ep1*ep2;         a[0][1] = ep1*f2;         a[0][2] = f1*ep2;         a[0][3] = f1*f2;
  a[1][0] = ep1*(1-ep2);     a[1][1] = ep1*(1-f2);     a[1][2] = (1-f1)*ep2;     a[1][3] = f1*(1-f2);
  a[2][0] = (1-ep1)*ep2;     a[2][1] = (1-ep1)*f2;     a[2][2] = (1-f1)*ep2;     a[2][3] = (1-f1)*f2;
  a[3][0] = (1-ep1)*(1-ep2); a[3][1] = (1-ep1)*(1-f2); a[3][2] = (1-f1)*(1-ep2); a[3][3] = (1-f1)*(1-f2);

  b[0] = Npp;                b[1] = Npf;               b[2] = Nfp;               b[3] = Nff;

  Jacobi(n,a,b,x, &count, &status);
  
  if(status == 2){
    cout << "no convergence in " << MAXIT << " iterations." << endl;
  }
  else{
    cout << "Solution in      " << count << " iterations:"<< endl;
    cout << "Nyy             =" << x[0] << endl;
    cout << "Nyj             =" << x[1] << endl;
    cout << "Njy             =" << x[2] << endl;
    cout << "Njj             =" << x[3] << endl;
    cout << "estimated purity:" << x[0]/(x[0]+x[1]+x[2]+x[3]) << endl;
  }
  

}
