#include <cstdio>
#include <cmath>
#include <vector>
#include <mwg/except.h>
#include "rktest.h"
#include "ksh/erk.h"

void test_erk1(){
  std::FILE* file;

  mwg_printd("euler");
  file = std::fopen("out/rk/rkeuler.txt","wb");
  test_method(file,kashiwa::rk16::euler_integrator());
  std::fclose(file);

  //
  // 2次公式
  //

  mwg_printd("mid");
  file = std::fopen("out/rk/rkmid.txt","wb");
  test_method(file,kashiwa::rk16::midpoint_integrator());
  std::fclose(file);

  mwg_printd("ralston");
  file = std::fopen("out/rk/rkral.txt","wb");
  test_method(file,kashiwa::rk16::ralston_integrator());
  std::fclose(file);

  mwg_printd("heun");
  file = std::fopen("out/rk/rkheun.txt","wb");
  test_method(file,kashiwa::rk16::heun_integrator());
  std::fclose(file);

  //
  // 3次公式
  //

  mwg_printd("runge3");
  file = std::fopen("out/rk/runge3.txt","wb");
  test_method(file,kashiwa::rk16::runge3_integrator());
  std::fclose(file);

  mwg_printd("heun3");
  file = std::fopen("out/rk/rkheun3.txt","wb");
  test_method(file,kashiwa::rk16::heun3_integrator());
  std::fclose(file);

  mwg_printd("ralston3");
  file = std::fopen("out/rk/rkral3.txt","wb");
  test_method(file,kashiwa::rk16::ralston3_integrator());
  std::fclose(file);

  mwg_printd("kutta3");
  file = std::fopen("out/rk/rkk3.txt","wb");
  test_method(file,kashiwa::rk16::kutta3_integrator());
  std::fclose(file);

  //
  // 4次公式
  //

  mwg_printd("gill");
  file = std::fopen("out/rk/rkgill.txt","wb");
  test_method(file,kashiwa::rk16::gill_integrator());
  std::fclose(file);

  mwg_printd("rk4");
  file = std::fopen("out/rk/rkrk4.txt","wb");
  test_method(file,kashiwa::rk16::rk4_integrator());
  std::fclose(file);

  mwg_printd("kutta38");
  file = std::fopen("out/rk/rkk38.txt","wb");
  test_method(file,kashiwa::rk16::kutta_3_8_integrator());
  std::fclose(file);

  //
  // 5次公式
  //

  mwg_printd("b5v1");
  file = std::fopen("out/rk/b5v1.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v1_integrator());
  std::fclose(file);

  mwg_printd("b5v2");
  file = std::fopen("out/rk/b5v2.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v2_integrator());
  std::fclose(file);

  mwg_printd("b5v3");
  file = std::fopen("out/rk/b5v3.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v3_integrator());
  std::fclose(file);

  //
  // 高次公式
  //

  mwg_printd("hammud6");
  file = std::fopen("out/rk/hammud6.txt","wb");
  test_method(file,kashiwa::rk16::hammud6_integrator());
  std::fclose(file);

  mwg_printd("shanks7");
  file = std::fopen("out/rk/shanks7.txt","wb");
  test_method(file,kashiwa::rk16::shanks7_integrator());
  std::fclose(file);

  mwg_printd("cv7");
  file = std::fopen("out/rk/rkcv7.txt","wb");
  test_method(file,kashiwa::rk16::cooper_verner7_integrator());
  std::fclose(file);

  mwg_printd("cv8");
  file = std::fopen("out/rk/rkcv8.txt","wb");
  test_method(file,kashiwa::rk16::cooper_verner8_integrator());
  std::fclose(file);
}
