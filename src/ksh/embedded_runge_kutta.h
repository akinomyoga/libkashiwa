// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_RK_EMBEDDED_RUNGE_KUTTA_H
#define KASHIWA_RK_EMBEDDED_RUNGE_KUTTA_H
#include <cstddef>
#include <cfloat>
#include "buffer.h"

namespace kashiwa{
namespace rk16{

  // DOP853: Dormand-Prince 8(5,3)
  //   from http://www.unige.ch/~hairer/prog/nonstiff/dop853.f (ハイラーの本の附録)
  struct dop853_integrator{
    static const int stage = 12;
    static const int order = 8;
    mutable working_buffer buffer;

    // buffer の状態
    //   before [  ? | k1 |  ? |  ? |  ? |  ? |  ? |  ? |  ? |  ? ]
    //   after  [  x | k1 | k5 | kC | x6 | x7 | x8 | k9 | kA | kB ]
    template<typename F>
    void _integrate8(
      double& time,double* __restrict__ value,std::size_t size,
      F const& f,double h,
      double atol,double rtol,double& _err,double& _stf
    ) const{
      double* __restrict__  x  = buffer.ptr();
      double* __restrict__  k1 = buffer.ptr()+size*1;
      double* __restrict__  k2 = buffer.ptr()+size*2;
      double* __restrict__& k3 = k2;
      double* __restrict__  k4 = buffer.ptr()+size*3;
      double* __restrict__& k5 = k3;
      double* __restrict__  k6 = buffer.ptr()+size*4;
      double* __restrict__  k7 = buffer.ptr()+size*5;
      double* __restrict__  k8 = buffer.ptr()+size*6;
      double* __restrict__  k9 = buffer.ptr()+size*7;
      double* __restrict__  kA = buffer.ptr()+size*8;
      double* __restrict__  kB = buffer.ptr()+size*9;
      double* __restrict__& kC = k4;

      static constexpr double c2 = 0.526001519587677318785587544488E-01;
      static constexpr double c3 = 0.789002279381515978178381316732E-01;
      static constexpr double c4 = 0.118350341907227396726757197510E+00;
      static constexpr double c5 = 0.281649658092772603273242802490E+00;
      static constexpr double c6 = 0.333333333333333333333333333333E+00;
      static constexpr double c7 = 0.25E+00;
      static constexpr double c8 = 0.307692307692307692307692307692E+00;
      static constexpr double c9 = 0.651282051282051282051282051282E+00;
      static constexpr double cA = 0.6E+00;
      static constexpr double cB = 0.857142857142857142857142857142E+00;
      static constexpr double cC = 1.0;

      static constexpr double a21 =  5.26001519587677318785587544488E-2;
      static constexpr double a31 =  1.97250569845378994544595329183E-2;
      static constexpr double a32 =  5.91751709536136983633785987549E-2;
      static constexpr double a41 =  2.95875854768068491816892993775E-2;
      static constexpr double a43 =  8.87627564304205475450678981324E-2;
      static constexpr double a51 =  2.41365134159266685502369798665E-1;
      static constexpr double a53 = -8.84549479328286085344864962717E-1;
      static constexpr double a54 =  9.24834003261792003115737966543E-1;
      static constexpr double a61 =  3.7037037037037037037037037037E-2;
      static constexpr double a64 =  1.70828608729473871279604482173E-1;
      static constexpr double a65 =  1.25467687566822425016691814123E-1;
      static constexpr double a71 =  3.7109375E-2;
      static constexpr double a74 =  1.70252211019544039314978060272E-1;
      static constexpr double a75 =  6.02165389804559606850219397283E-2;
      static constexpr double a76 = -1.7578125E-2;
      static constexpr double a81 =  3.70920001185047927108779319836E-2;
      static constexpr double a84 =  1.70383925712239993810214054705E-1;
      static constexpr double a85 =  1.07262030446373284651809199168E-1;
      static constexpr double a86 = -1.53194377486244017527936158236E-2;
      static constexpr double a87 =  8.27378916381402288758473766002E-3;
      static constexpr double a91 =  6.24110958716075717114429577812E-1;
      static constexpr double a94 = -3.36089262944694129406857109825;
      static constexpr double a95 = -8.68219346841726006818189891453E-1;
      static constexpr double a96 =  2.75920996994467083049415600797E+1;
      static constexpr double a97 =  2.01540675504778934086186788979E+1;
      static constexpr double a98 = -4.34898841810699588477366255144E+1;
      static constexpr double aA1 =  4.77662536438264365890433908527E-1;
      static constexpr double aA4 = -2.48811461997166764192642586468;
      static constexpr double aA5 = -5.90290826836842996371446475743E-1;
      static constexpr double aA6 =  2.12300514481811942347288949897E+1;
      static constexpr double aA7 =  1.52792336328824235832596922938E+1;
      static constexpr double aA8 = -3.32882109689848629194453265587E+1;
      static constexpr double aA9 = -2.03312017085086261358222928593E-2;
      static constexpr double aB1 = -9.3714243008598732571704021658E-1;
      static constexpr double aB4 =  5.18637242884406370830023853209;
      static constexpr double aB5 =  1.09143734899672957818500254654;
      static constexpr double aB6 = -8.14978701074692612513997267357;
      static constexpr double aB7 = -1.85200656599969598641566180701E+1;
      static constexpr double aB8 =  2.27394870993505042818970056734E+1;
      static constexpr double aB9 =  2.49360555267965238987089396762;
      static constexpr double aBA = -3.0467644718982195003823669022;
      static constexpr double aC1 =  2.27331014751653820792359768449;
      static constexpr double aC4 = -1.05344954667372501984066689879E+1;
      static constexpr double aC5 = -2.00087205822486249909675718444;
      static constexpr double aC6 = -1.79589318631187989172765950534E+1;
      static constexpr double aC7 =  2.79488845294199600508499808837E+1;
      static constexpr double aC8 = -2.85899827713502369474065508674;
      static constexpr double aC9 = -8.87285693353062954433549289258;
      static constexpr double aCA =  1.23605671757943030647266201528E+1;
      static constexpr double aCB =  6.43392746015763530355970484046E-1;

      static constexpr double b1 =  5.42937341165687622380535766363E-2;
      static constexpr double b6 =  4.45031289275240888144113950566;
      static constexpr double b7 =  1.89151789931450038304281599044;
      static constexpr double b8 = -5.8012039600105847814672114227;
      static constexpr double b9 =  3.1116436695781989440891606237E-1;
      static constexpr double bA = -1.52160949662516078556178806805E-1;
      static constexpr double bB =  2.01365400804030348374776537501E-1;
      static constexpr double bC =  4.47106157277725905176885569043E-2;

      static constexpr double er1 =  0.1312004499419488073250102996E-01;
      static constexpr double er6 = -0.1225156446376204440720569753E+01;
      static constexpr double er7 = -0.4957589496572501915214079952E+00;
      static constexpr double er8 =  0.1664377182454986536961530415E+01;
      static constexpr double er9 = -0.3503288487499736816886487290E+00;
      static constexpr double erA =  0.3341791187130174790297318841E+00;
      static constexpr double erB =  0.8192320648511571246570742613E-01;
      static constexpr double erC = -0.2235530786388629525884427845E-01;

      static constexpr double bhh1 = 0.244094488188976377952755905512E+00;
      static constexpr double bhh9 = 0.733846688281611857341361741547E+00;
      static constexpr double bhhC = 0.220588235294117647058823529412E-01;

      // k1: この関数を呼び出す前に設定しておく (FSAL)

      // k2
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a21*k1[i]);
      f(k2,time+c2*h,x);

      // k3 := k2
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a31*k1[i]+a32*k2[i]);
      f(k3,time+c3*h,x);

      // k4
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a41*k1[i]+a43*k3[i]);
      f(k4,time+c4*h,x);

      // k5 := k3
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a51*k1[i]+a53*k3[i]+a54*k4[i]);
      f(k5,time+c5*h,x);

      // k6
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a61*k1[i]+a64*k4[i]+a65*k5[i]);
      f(k6,time+c6*h,x);

      // k7
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a71*k1[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
      f(k7,time+c7*h,x);

      // k8
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a81*k1[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
      f(k8,time+c8*h,x);

      // k9
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(a91*k1[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+a98*k8[i]);
      f(k9,time+c9*h,x);

      // kA
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(aA1*k1[i]+aA4*k4[i]+aA5*k5[i]+aA6*k6[i]+aA7*k7[i]+aA8*k8[i]+aA9*k9[i]);
      f(kA,time+cA*h,x);

      // kB
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(aB1*k1[i]+aB4*k4[i]+aB5*k5[i]+aB6*k6[i]+aB7*k7[i]+aB8*k8[i]+aB9*k9[i]+aBA*kA[i]);
      f(kB,time+cB*h,x);

      // kC
      for(std::size_t i=0;i<size;i++)
        x[i] = value[i] + h*(aC1*k1[i]+aC4*k4[i]+aC5*k5[i]+aC6*k6[i]+aC7*k7[i]+aC8*k8[i]+aC9*k9[i]+aCA*kA[i]+aCB*kB[i]);
      f(kC,time+cC*h,x);

      // increment
      double err1 = 0.0, err2 = 0.0, err3 = 0.0;
      for(std::size_t i=0;i<size;i++){
        double const slope = b1*k1[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+bA*kA[i]+bB*kB[i]+bC*kC[i];
        double const _y = value[i] + h*slope;
        double const sk = atol+rtol*std::max(value[i],_y);
        double const e1 = er1*k1[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+erA*kA[i]+erB*kB[i]+erC*kC[i];
        double const e2 = slope-bhh1*k1[i]-bhh9*k9[i]-bhhC*kC[i];

        x[i] = _y;
        err1 += (e1/sk)*(e1/sk);
        err2 += (e2/sk)*(e2/sk);
        err3 += (_y-x[i])*(_y-x[i]);
      }

      double deno = err1+0.01*err2;
      if(deno<=0.0) deno = 1.0;
      _err = std::abs(h)*err1/std::sqrt(size*deno);
      _stf = err3;
    }

    struct stat_t{
      int nfcn   { 0 };
      int nstep  { 0 };
      int naccpt { 0 };
      int nrejct { 0 };
    };

    struct param_t{
      double atol {1e-13};
      double rtol {1e-13};
      double beta {0.0};
      double fac1 {0.0};
      double fac2 {0.0};
      double safe {0.0};
      double hmax {0.0};
      double step {0.0};
      std::ptrdiff_t nmax {100000};
    };

    // buffer の状態
    //   before [  x | k1 | kD | kC | k6 | k7 | k8 | k9 | kA | kB ]
    //   after  [  x | k1 | kD | kC | xE | xF | xG | kE | kF | kG ]
    template<typename F>
    void _dense_output_initialize(
      double& time,double const* __restrict__ value,std::size_t size,
      F const& f,double h,stat_t& stat,int* icomp,std::size_t nrd,double* __restrict__ cont,double& hout
    ) const{
      double* __restrict__  x = buffer.ptr();
      double* __restrict__  k1 = buffer.ptr()+size*1;
      double* __restrict__  k6 = buffer.ptr()+size*4;
      double* __restrict__  k7 = buffer.ptr()+size*5;
      double* __restrict__  k8 = buffer.ptr()+size*6;
      double* __restrict__  k9 = buffer.ptr()+size*7;
      double* __restrict__  kA = buffer.ptr()+size*8;
      double* __restrict__  kB = buffer.ptr()+size*9;
      double* __restrict__  kC = buffer.ptr()+size*3;
      double* __restrict__  kD = buffer.ptr()+size*2;

      double* __restrict__& xE = k6;
      double* __restrict__& xF = k7;
      double* __restrict__& xG = k8;
      double* __restrict__& kE = k9;
      double* __restrict__& kF = kA;
      double* __restrict__& kG = kB;

      static constexpr double cE  = 0.1E+00;
      static constexpr double cF  = 0.2E+00;
      static constexpr double cG  = 0.777777777777777777777777777778E+00;

      static constexpr double aE1 =  5.61675022830479523392909219681E-2;
      static constexpr double aE7 =  2.53500210216624811088794765333E-1;
      static constexpr double aE8 = -2.46239037470802489917441475441E-1;
      static constexpr double aE9 = -1.24191423263816360469010140626E-1;
      static constexpr double aEA =  1.5329179827876569731206322685E-1;
      static constexpr double aEB =  8.20105229563468988491666602057E-3;
      static constexpr double aEC =  7.56789766054569976138603589584E-3;
      static constexpr double aED = -8.298E-3;
      static constexpr double aF1 =  3.18346481635021405060768473261E-2;
      static constexpr double aF6 =  2.83009096723667755288322961402E-2;
      static constexpr double aF7 =  5.35419883074385676223797384372E-2;
      static constexpr double aF8 = -5.49237485713909884646569340306E-2;
      static constexpr double aFB = -1.08347328697249322858509316994E-4;
      static constexpr double aFC =  3.82571090835658412954920192323E-4;
      static constexpr double aFD = -3.40465008687404560802977114492E-4;
      static constexpr double aFE =  1.41312443674632500278074618366E-1;
      static constexpr double aG1 = -4.28896301583791923408573538692E-1;
      static constexpr double aG6 = -4.69762141536116384314449447206;
      static constexpr double aG7 =  7.68342119606259904184240953878;
      static constexpr double aG8 =  4.06898981839711007970213554331;
      static constexpr double aG9 =  3.56727187455281109270669543021E-1;
      static constexpr double aGD = -1.39902416515901462129418009734E-3;
      static constexpr double aGE =  2.9475147891527723389556272149;
      static constexpr double aGF = -9.15095847217987001081870187138;

      static constexpr double d41 = -0.84289382761090128651353491142E+01;
      static constexpr double d46 =  0.56671495351937776962531783590E+00;
      static constexpr double d47 = -0.30689499459498916912797304727E+01;
      static constexpr double d48 =  0.23846676565120698287728149680E+01;
      static constexpr double d49 =  0.21170345824450282767155149946E+01;
      static constexpr double d4A = -0.87139158377797299206789907490E+00;
      static constexpr double d4B =  0.22404374302607882758541771650E+01;
      static constexpr double d4C =  0.63157877876946881815570249290E+00;
      static constexpr double d4D = -0.88990336451333310820698117400E-01;
      static constexpr double d4E =  0.18148505520854727256656404962E+02;
      static constexpr double d4F = -0.91946323924783554000451984436E+01;
      static constexpr double d4G = -0.44360363875948939664310572000E+01;

      static constexpr double d51 =  0.10427508642579134603413151009E+02;
      static constexpr double d56 =  0.24228349177525818288430175319E+03;
      static constexpr double d57 =  0.16520045171727028198505394887E+03;
      static constexpr double d58 = -0.37454675472269020279518312152E+03;
      static constexpr double d59 = -0.22113666853125306036270938578E+02;
      static constexpr double d5A =  0.77334326684722638389603898808E+01;
      static constexpr double d5B = -0.30674084731089398182061213626E+02;
      static constexpr double d5C = -0.93321305264302278729567221706E+01;
      static constexpr double d5D =  0.15697238121770843886131091075E+02;
      static constexpr double d5E = -0.31139403219565177677282850411E+02;
      static constexpr double d5F = -0.93529243588444783865713862664E+01;
      static constexpr double d5G =  0.35816841486394083752465898540E+02;

      static constexpr double d61 =  0.19985053242002433820987653617E+02;
      static constexpr double d66 = -0.38703730874935176555105901742E+03;
      static constexpr double d67 = -0.18917813819516756882830838328E+03;
      static constexpr double d68 =  0.52780815920542364900561016686E+03;
      static constexpr double d69 = -0.11573902539959630126141871134E+02;
      static constexpr double d6A =  0.68812326946963000169666922661E+01;
      static constexpr double d6B = -0.10006050966910838403183860980E+01;
      static constexpr double d6C =  0.77771377980534432092869265740E+00;
      static constexpr double d6D = -0.27782057523535084065932004339E+01;
      static constexpr double d6E = -0.60196695231264120758267380846E+02;
      static constexpr double d6F =  0.84320405506677161018159903784E+02;
      static constexpr double d6G =  0.11992291136182789328035130030E+02;

      static constexpr double d71 = -0.25693933462703749003312586129E+02;
      static constexpr double d76 = -0.15418974869023643374053993627E+03;
      static constexpr double d77 = -0.23152937917604549567536039109E+03;
      static constexpr double d78 =  0.35763911791061412378285349910E+03;
      static constexpr double d79 =  0.93405324183624310003907691704E+02;
      static constexpr double d7A = -0.37458323136451633156875139351E+02;
      static constexpr double d7B =  0.10409964950896230045147246184E+03;
      static constexpr double d7C =  0.29840293426660503123344363579E+02;
      static constexpr double d7D = -0.43533456590011143754432175058E+02;
      static constexpr double d7E =  0.96324553959188282948394950600E+02;
      static constexpr double d7F = -0.39177261675615439165231486172E+02;
      static constexpr double d7G = -0.14972683625798562581422125276E+03;

      for(int j=0; j<nrd; j++){
        int const i=icomp[j]; // 出力する変数の index
        cont[j]=value[i];
        double const ydiff=x[i]-value[i];
        cont[j+nrd]=ydiff;
        double const bspl=h*k1[i]-ydiff;
        cont[j+nrd*2]=bspl;
        cont[j+nrd*3]=ydiff-h*kD[i]-bspl;
        cont[j+nrd*4]=d41*k1[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+d49*k9[i]+d4A*kA[i]+d4B*kB[i]+d4C*kC[i];
        cont[j+nrd*5]=d51*k1[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+d59*k9[i]+d5A*kA[i]+d5B*kB[i]+d5C*kC[i];
        cont[j+nrd*6]=d61*k1[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+d69*k9[i]+d6A*kA[i]+d6B*kB[i]+d6C*kC[i];
        cont[j+nrd*7]=d71*k1[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+d79*k9[i]+d7A*kA[i]+d7B*kB[i]+d7C*kC[i];
      }

      for(std::size_t i=0;i<size;i++){
        double const _xE = value[i] + h*(aE1*k1[i]+aE7*k7[i]+aE8*k8[i]+aE9*k9[i]+aEA*kA[i]+aEB*kB[i]+aEC*kC[i]+aED*kD[i]);
        double const _xF = value[i] + h*(aF1*k1[i]+aF6*k6[i]+aF7*k7[i]+aF8*k8[i]+aFB*kB[i]+aFC*kC[i]+aFD*kD[i]);
        double const _xG = value[i] + h*(aG1*k1[i]+aG6*k6[i]+aG7*k7[i]+aG8*k8[i]+aG9*k9[i]+aGD*kD[i]);
        /* xE == k6 */ k6[i] = _xE;
        /* xF == k7 */ k7[i] = _xF;
        /* xG == k8 */ k8[i] = _xG;
      }

      f(kE,time+cE*h,xE);

      for(std::size_t i=0;i<size;i++) xF[i] += h*aFE*kE[i];
      f(kF,time+cF*h,xF);

      for(std::size_t i=0;i<size;i++) xG[i] += h*(aGE*kE[i]+aGF*kF[i]);
      f(kG,time+cG*h,xG);

      stat.nfcn+=3;
      for(int j=0; j<nrd; j++){
        int const i=icomp[j];
        cont[j+nrd*4]=h*(cont[j+nrd*4]+d4D*kD[i]+d4E*kE[i]+d4F*kF[i]+d4G*kG[i]);
        cont[j+nrd*5]=h*(cont[j+nrd*5]+d5D*kD[i]+d5E*kE[i]+d5F*kF[i]+d5G*kG[i]);
        cont[j+nrd*6]=h*(cont[j+nrd*6]+d6D*kD[i]+d6E*kE[i]+d6F*kF[i]+d6G*kG[i]);
        cont[j+nrd*7]=h*(cont[j+nrd*7]+d7D*kD[i]+d7E*kE[i]+d7F*kF[i]+d7G*kG[i]);
      }
      hout=h;
    }

    template<typename F>
    double _determine_initial_step(
      double time,double* __restrict__ value,std::size_t size,F const& f,
      int bwd,double atol,double rtol,double hmax
    ) const{
      double* __restrict__ const x  = buffer.ptr();
      double* __restrict__ const k1 = buffer.ptr()+size*1;
      double* __restrict__ const k2 = buffer.ptr()+size*2;

      // compute a first guess for explicit euler as:
      //   h1 = 0.01 * NORM(value) / NORM(k1).
      double dnf=0.0, dny=0.0;
      for(std::size_t i = 0; i<size; i++){
        double const sk=atol+rtol*std::abs(value[i]);
        dnf+=(k1[i]/sk)*(k1[i]/sk);
        dny+=(value[i]/sk)*(value[i]/sk);
      }

      double h1;
      if(dnf<=1e-10||dny<=1e-10)
        h1=1.0e-6;
      else
        h1=std::sqrt(dny/dnf)*0.01;
      h1=std::min(h1,hmax);

      for(std::size_t i = 0; i<size; i++)
        x[i]=value[i]+bwd*h1*k1[i];
      f(k2,time+bwd*h1,x);

      // 二階微分ノルム |f'|
      double norm2=0.0;
      for(std::size_t i = 0; i<size; i++){
        double const sk=atol+rtol*abs(value[i]);
        double const z =(k2[i]-k1[i])/sk;
        norm2+=z*z;
      }
      norm2=std::sqrt(norm2)/h1;

      // step size `h2' is computed such that:
      //   h2^order * max{norm(f), norm(f')} = 0.01
      double const der12=std::max(norm2,std::sqrt(dnf));
      double h2;
      if(der12<=1e-15)
        h2=std::max(1e-6,std::abs(h1)*1e-3);
      else
        h2=std::pow(0.01/der12, 1.0/order);

      return bwd*std::min(std::min(100*h1,h2),hmax);
    }

  public:
    template<typename F>
    void operator()(double& time,double* __restrict__ value,std::size_t size,F const& f,double h) const{
      buffer.ensure(10*size);
      double const atol = 1e-12;
      double const rtol = 1e-12;
      double err,stf;

      double* __restrict__  x = buffer.ptr();
      double* __restrict__  k1 = buffer.ptr()+size*1;
      f(k1,time,value);
      this->_integrate8(
        time,value,size,f,h,
        atol,rtol,err,stf
      );

      for(std::size_t i=0;i<size;i++)
        value[i] = x[i];
      time+=h;
    }

    template<typename F>
    void integrate(
      double& time,double* __restrict__ value,std::size_t size,F const& f,
      double timeN,stat_t& stat,param_t const& params
    ) const{
      buffer.ensure(10*size);
      double const beta  = kashiwa::clamp(params.beta, 0.0, 0.2); // e.g. 0.04
      double const safe  = params.safe==0.0?0.9: kashiwa::clamp(params.safe, 1e-4, 1.0);
      double const facc1 = 1.0/(params.fac1==0.0?0.333: params.fac1);
      double const facc2 = 1.0/(params.fac2==0.0?6.000: params.fac2);
      double const expo1 = 1.0/8.0-beta*0.2;
      double const hmax  = std::abs(params.hmax==0.0?timeN-time: params.hmax);
      int    const bwd   = time<timeN?1: -1;
      double const rtol  = std::abs(params.rtol);
      double const atol  = std::abs(params.atol);
      std::ptrdiff_t const nmax = params.nmax;

      double* __restrict__ const x  = buffer.ptr();
      double* __restrict__ const k1 = buffer.ptr()+size*1;
      double* __restrict__ const kD = buffer.ptr()+size*2;

      f(k1,time,value);
      stat.nfcn++;

      double h = bwd*std::abs(params.step);
      if(h == 0.0){
        this->_determine_initial_step(time,value,size,f,bwd,atol,rtol,hmax);
        stat.nfcn++;
      }

      double facold=1e-4;
      bool reject=false, last=false;
      // double hlamb=0.0;
      // int iasti=0;
      for(;;stat.nstep++){
        mwg_check(nmax<0||stat.nstep<nmax,"収束しません time = %g, h = %g at step#%d",time,h,stat.nstep);
        mwg_check(0.1*std::abs(h)>std::abs(time)*DBL_EPSILON,"時刻桁落ち time = %g, h = %g",time,h);

        if((time+1.01*h-timeN)*bwd>0.0){
          h=timeN-time;
          last=true;
        }

        //mwg_printd("time = %g, h = %g",time,h);
        double err,stf;
        this->_integrate8(
          time,value,size,f,h,
          atol,rtol,err,stf
        );
        stat.nfcn+=11;

        double const fac11=std::pow(err,expo1);
        double const fac = kashiwa::clamp(fac11/(std::pow(facold,beta)*safe), facc2, facc1);
        double hnew = h/fac;
        if(err > 1.0){
          h /= std::min(facc1, fac11/safe);
          reject=true;
          last=false;
          if(stat.naccpt>=1) stat.nrejct++;
        }else{
          stat.naccpt++;
          facold=std::max(err,1e-4);
          f(kD,time+h,x);
          stat.nfcn++;

          // stiffness detection
          // if(stat.naccpt%nstiff==0||iasti>0){
          //   if(stf>0.0){
          //     double stnum=0.0;
          //     for(std::size_t i=0;i<size;i++)
          //       stnum=stnum+(kD[i]-kC[i])**2;
          //     hlamb=std::abs(h)*std::sqrt(stnum/stf);
          //   }
          //   if(hlamb>6.1){
          //     nonsti=0;
          //     iasti++;
          //     if(iasti==15)
          //       mwg_printd("the problem seems to become stiff at time = %g",time);
          //   }else{
          //     nonsti++;
          //     if(nonsti==6) iasti=0;
          //   }
          // }

          // 密出力

          // if(iout==2 || event)
          //   _dense_output_initialize(...);

          // if(iout==1||iout==2||event){
          //   solout(naccpt+1,time,time+h,value,size,cont,icomp,nrd,irtrn,xout);
          //   if (irtrn<0) goto 79; // error report
          // }

          for(std::size_t i=0;i<size;i++){
            k1[i] = kD[i];
            value[i] = x[i];
          }

          time+=h;
          if(last)return;

          if(std::abs(hnew)>hmax) hnew=bwd*hmax;
          if(reject) hnew=bwd*std::min(std::abs(hnew),std::abs(h));
          reject=false;

          h=hnew;
        }
      }
    }
  };

}
}
#endif
