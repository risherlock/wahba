#include "esoq2.h"
#include <stdio.h>

int main()
{
  const int n = 2;
  esoq2_vec_t vb[2], vr[2];

  /* Reference vectors */

  // Magnetometer
  vr[0].x = 0.259467526530057;
  vr[0].y = 0.664468321860749;
  vr[0].z = 0.700826977163361;

  // Accelerometer
  vr[1].x = 0.492829407411352;
  vr[1].y = 0.610522688331208;
  vr[1].z = 0.619984856446841;

  /* Measurements */

  // Magnetometer
  vb[0].x = 3.105722495759851e-03;
  vb[0].y = 9.827958004270893e-01;
  vb[0].z = 2.014309934927343e-01;

  // Accelerometer
  vb[1].x = 2.495092383950536e-01;
  vb[1].y = 9.391322928827616e-01;
  vb[1].z = 2.402395094056679e-01;

  double w[] = {0.01, 0.01};
  quaternion_t q = esoq2_update(vb, vr, w, n);

  printf("q0: %f\n", q.q0);
  printf("q1: %f\n", q.q1);
  printf("q2: %f\n", q.q2);
  printf("q3: %f\n", q.q3);

  return 0;
}

// q_hat
//    3.112751477015743e-01
//    1.717597371574912e-01
//   -5.136381752360104e-02
//    9.332567349686552e-01
