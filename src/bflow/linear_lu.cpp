#include <cstdio>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
//#include <mwg/except.h>
#include "linear_lu.h"

namespace{

  class LUDecomposer{
    int N;
    double* arr;
    int* imap;

    double& A(int i,int j){
      return arr[imap[i]*N+j];
    }

    double determine_pivot(int i){
      int imax=i;
      double vmax=std::abs(A(i,i));
      for(int icand=i+1;icand<N;icand++){
        double const vcand=std::abs(A(icand,i));
        if(vcand>vmax){
          imax=icand;
          vmax=vcand;
        }
      }
      if(imax>i){
        using namespace std;
        swap(imap[i],imap[imax]);
      }
      //mwg_assert_nothrow(vmax!=0);
      return A(i,i); // pivot 交換後の対角成分
    }

  public:
    void decompose(double* arr,int* imap){
      this->arr=arr;
      this->imap=imap;

      for(int i=0;i<N;i++)imap[i]=i;

      for(int i=0;i<N;i++){
        double const scal=1.0/this->determine_pivot(i);

        for(int ii=i+1;ii<N;ii++){
          double const l=A(ii,i)*=scal;
          for(int jj=i+1;jj<N;jj++){
            double const u=A(i,jj);
            A(ii,jj)-=l*u;
          }
        }
      }
    }

  public:
    LUDecomposer(int N):N(N){}
  };

  class LUEquationSolver{
    int N;
    double* arr;
    double* vec;
    int* imap;
    std::vector<double> vtmp;

    double& A(int i,int j){
      return arr[imap[i]*N+j];
    }
    double& b(int i){
      return vec[imap[i]];
    }

    void forward_substitute(){
      for(int i=0;i<N;i++){
        double& vv(vtmp[i]=b(i));
        for(int j=0;j<i;j++)
          vv-=A(i,j)*vtmp[j];
      }
    }
    void backward_substitute(){
      for(int i=N-1;i>=0;i--){
        double& vv(vec[i]=vtmp[i]);
        for(int j=i+1;j<N;j++)
          vv-=A(i,j)*vec[j];
        vv/=A(i,i);
      }
    }

  public:
    void solve(double const* arr,int const* imap,double* vec){
      this->arr=const_cast<double*>(arr);
      this->imap=const_cast<int*>(imap);
      this->vec=vec;
      this->forward_substitute();
      this->backward_substitute();
    }

  public:
    LUEquationSolver(int N):N(N){
      this->vtmp.resize(N,0.0);
    }
  };


  // template<int N>
  // void SolveLinearEquationLU(double* arr,double* vec){
  //   // 取り敢えずの実装:
  //   //   arr を破壊的に使用する
  //   //   vec ベクトルを指定し結果を格納する。
  //   LinearEquationSolver_LU(N,arr,vec);
  // }

}

namespace idt{
namespace rfh{
  void LUDecompose(int N,double* arr,int* imap){
    LUDecomposer calculator(N);
    calculator.decompose(arr,imap);
  }
  void LUEquationSolve(int N,double const* arr,int const* imap,double* vec){
    LUEquationSolver calculator(N);
    calculator.solve(arr,imap,vec);
  }
  void SolveLinearEquationLU(std::size_t N,double* arr,double* vec){
    std::vector<int> imap((std::size_t)N);
    LUDecompose(N,arr,&imap[0]);
    LUEquationSolve(N,arr,&imap[0],vec);
  }
  void SolveLinearEquation(int N,double* arr,double* vec){
    SolveLinearEquationLU(N,arr,vec);
  }

}
}

// int main(){
//   return 0;
// }
