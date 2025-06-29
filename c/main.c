#include "esoq2.h"
#include <stdio.h>
#include <math.h>

int main()
{
  const int n = 5;          // Number of sensor pairs
  double w[n];              // Sensor weights
  esoq2_vec_t vr[n], vb[n]; // Reference and measurement vectors

  // Reference vectors
  double r_data[3][5] = {
      {6.986305305048354e-02, 2.401869784514148e-01, 8.391504991198349e-01, 3.600210349220398e-01, 7.027533847857507e-01},
      {6.247648709321313e-01, 8.911751620367423e-01, 4.756028241780341e-01, 6.504202271422250e-01, 7.105712111509511e-01},
      {7.776811749474347e-01, 3.848597743999308e-01, 2.638719262460676e-01, 6.688335985414607e-01, 3.501762492861636e-02}};

  // Measurement vectors in body frame
  double b_data[3][5] = {
      {3.459504954471201e-01, 6.082988033461165e-01, 9.740580971739062e-01, 6.206425793515392e-01, 9.401449260928908e-01},
      {4.394569564964670e-01, 6.582609813328052e-01, 5.855150106905321e-02, 3.632267836446652e-01, 3.581766519779687e-01},
      {8.349960702405268e-01, 4.629623909219600e-01, 2.537410682802650e-01, 7.110326219020270e-01, 5.938731418624136e-02}};

  // Initialize vectors
  for (int i = 0; i < n; i++)
  {
    vr[i].x = r_data[0][i];
    vr[i].y = r_data[1][i];
    vr[i].z = r_data[2][i];

    vb[i].x = b_data[0][i];
    vb[i].y = b_data[1][i];
    vb[i].z = b_data[2][i];

    w[i] = 1.0;
  }

  quaternion_t q_hat = esoq2_update(vb, vr, w, n);
  quaternion_t q_truth = {-5.967649364834088e-02, -2.564977525174380e-02, 2.116648784998131e-01, 9.751814109923517e-01};

  /* Results */
  printf("Estimated quaternion:\n");
  printf("q0: %.15f\n", q_hat.q0);
  printf("q1: %.15f\n", q_hat.q1);
  printf("q2: %.15f\n", q_hat.q2);
  printf("q3: %.15f\n\n", q_hat.q3);

  printf("True quaternion:\n");
  printf("q0: %.15f\n", q_truth.q3);
  printf("q1: %.15f\n", q_truth.q0);
  printf("q2: %.15f\n", q_truth.q1);
  printf("q3: %.15f\n", q_truth.q2);

  return 0;
}

// Estimated quaternion:
// q0: 0.975450636761736
// q1: -0.060016533883873
// q2: -0.026049455305635
// q3: 0.210274812519163

// True quaternion:
// q0: -0.059676493648341
// q1: -0.025649775251744
// q2: 0.211664878499813
// q3: 0.975181410992352
