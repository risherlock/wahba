#include <math.h>

#include "esoq2.h"

static inline double det_3x3(double m[3][3])
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

// Compute greatest eigenvalue of traceless symmetric 4x4 matrix
double get_lambda_max(const double K[4][4], const int n)
{
  // Cache unique elements of the symmetric matrix
  const double K00 = K[0][0], K01 = K[0][1], K02 = K[0][2], K03 = K[0][3];
  const double K11 = K[1][1], K12 = K[1][2], K13 = K[1][3];
  const double K22 = K[2][2], K23 = K[2][3];
  const double K33 = K[3][3];

  // Cache frequently used square terms
  const double K03_sq = K03 * K03;
  const double K13_sq = K13 * K13;
  const double K23_sq = K23 * K23;
  const double K12_sq = K12 * K12;
  const double K01_sq = K01 * K01;
  const double K02_sq = K02 * K02;

  // Trace and adjugate terms. Eqn.(7)
  const double trB = K33;
  const double tradB = // trace(adj(B + tr(B)))
      (K11 + trB) * (K22 + trB) - K12_sq +
      (K00 + trB) * (K22 + trB) - K02_sq +
      (K00 + trB) * (K11 + trB) - K01_sq;
  const double b = -2.0f * trB * trB + tradB - (K03_sq + K13_sq + K23_sq);

  // Determinant of symmetric 4x4 matrix
  const double d =
      K03_sq * K12_sq - K00 * K33 * K12_sq - 2.0 * K01 * K03 * K23 * K12 +
      2.0 * K01 * K02 * K33 * K12 + 2.0 * K00 * K23 * K12 * K13 -
      2.0 * K02 * K03 * K12 * K13 - K11 * K00 * K23_sq + K01_sq * K23_sq +
      2.0 * K11 * K02 * K03 * K23 - K11 * K02_sq * K33 + K02_sq * K13_sq -
      2.0 * K01 * K02 * K23 * K13 + K11 * K00 * K33 * K22 - K00 * K13_sq * K22 -
      K01_sq * K33 * K22 + 2.0 * K01 * K03 * K13 * K22 - K11 * K03_sq * K22;

  // Two quaternions
  if (n == 2)
  {
    const double sqrt_d = sqrt(d);
    const double g3 = sqrt(2 * sqrt_d - b);  // Eqn.(15)
    const double g4 = sqrt(-2 * sqrt_d - b); // Eqn.(15)
    return (g3 + g4) / 2.0f;                 // Eqn.(14)
  }

  // trace(adj(K))
  double m1[3][3] = {{K11, K12, K13}, {K12, K22, K23}, {K13, K23, K33}};
  double m2[3][3] = {{K00, K02, K03}, {K02, K22, K23}, {K03, K23, K33}};
  double m3[3][3] = {{K00, K01, K03}, {K01, K11, K13}, {K03, K13, K33}};
  double m4[3][3] = {{K00, K01, K02}, {K01, K11, K12}, {K02, K12, K22}};
  const double tradK = det_3x3(m1) + det_3x3(m2) + det_3x3(m3) + det_3x3(m4);
  const double c = -tradK; // Eqn.(7)

  // Eqn.(10)
  const double p = b * b / 9.0 + 4.0 * d / 3.0;
  const double q = b * b * b / 27.0 - 4.0 * d * b / 3.0 + c * c / 2.0;

  // Eigenvalue solution
  const double arg = q / pow(p, 1.5);
  const double theta = acos(fmax(-1.0, fmin(1.0, arg)));
  const double u1 = 2.0 * sqrt(p) * cos(theta / 3.0) + b / 3.0; // Eqn.(9)
  const double g1 = sqrt(u1 - b);                               // Eqn.(12)
  const double g2 = -2.0 * sqrt(u1 * u1 - 4.0 * d);             // Eqn.(12)
  return 0.5f * (g1 + sqrt(-u1 - b - g2));                      // Eqn.(13)
}

quaternion_t esoq2_update(const esoq2_vec_t *vb, const esoq2_vec_t *vi, double *w, const int n)
{
  double B[3][3] = {{0.0}, {0.0}, {0.0}};

  for (int i = 0; i < n; i++)
  {
    B[0][0] += w[i] * vb[i].x * vi[i].x;
    B[0][1] += w[i] * vb[i].x * vi[i].y;
    B[0][2] += w[i] * vb[i].x * vi[i].z;

    B[1][0] += w[i] * vb[i].y * vi[i].x;
    B[1][1] += w[i] * vb[i].y * vi[i].y;
    B[1][2] += w[i] * vb[i].y * vi[i].z;

    B[2][0] += w[i] * vb[i].z * vi[i].x;
    B[2][1] += w[i] * vb[i].z * vi[i].y;
    B[2][2] += w[i] * vb[i].z * vi[i].z;
  }

  const double trB = B[0][0] + B[1][1] + B[2][2];
  const double z[3] = {B[1][2] - B[2][1], B[2][0] - B[0][2], B[0][1] - B[1][0]};

  // S = B + B' - I * trB
  double S00 = 2 * B[0][0] - trB;
  double S11 = 2 * B[1][1] - trB;
  double S22 = 2 * B[2][2] - trB;
  const double S01 = B[0][1] + B[1][0];
  const double S02 = B[0][2] + B[2][0];
  const double S12 = B[1][2] + B[2][1];

  // Compute largest eigenvalue of K
  const double K[4][4] = {{S00, S01, S02, z[0]}, {S01, S11, S12, z[1]}, {S02, S12, S22, z[2]}, {z[0], z[1], z[2], trB}};
  const double lambda = get_lambda_max(K, n);

  // Eqn.(18)
  S00 = S00 - lambda;
  S11 = S11 - lambda;
  S22 = S22 - lambda;
  const double t = trB - lambda;

  // Entries of matrix M in Eqn.(20) computed using Eqn.(19)
  const double ma = S00 * t - z[0] * z[0];
  const double mb = S11 * t - z[1] * z[1];
  const double mc = S22 * t - z[2] * z[2];
  const double mx = S01 * t - z[0] * z[1];
  const double my = S02 * t - z[0] * z[2];
  const double mz = S12 * t - z[1] * z[2];

  double e[3] = {mb * mc - mz * mz, ma * mc - my * my, ma * mb - mx * mx};
  double e_abs[3] = {fabs(e[0]), fabs(e[1]), fabs(e[2])};

  // Choose optimal principal axis from three choices. Eqn.(21)
  if ((e_abs[0] > e_abs[1]) && (e_abs[0] > e_abs[2]))
  {
    e[1] = my * mz - mx * mc;
    e[2] = mx * mz - my * mb;
  }
  else if (e_abs[1] > e_abs[2])
  {
    e[0] = my * mz - mx * mc;
    e[2] = mx * my - mz * ma;
  }
  else
  {
    e[0] = mx * mz - my * mb;
    e[1] = mx * my - mz * ma;
  }

  // Normalize e
  const double norm_e_inv = 1.0 / sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  e[0] *= norm_e_inv;
  e[1] *= norm_e_inv;
  e[2] *= norm_e_inv;

  // Find maximum absolute value among the elements of z and t
  const double abs_z0 = fabs(z[0]);
  const double abs_z1 = fabs(z[1]);
  const double abs_z2 = fabs(z[2]);
  const double abs_t = fabs(t);

  double xk = 0.0;
  double yk = 0.0;

  if ((abs_z0 > abs_z1) && (abs_z0 > abs_z2) && (abs_z0 > abs_t))
  {
    xk = z[0];
    yk = S00 * e[0] + S01 * e[1] + S02 * e[2];
  }
  else if ((abs_z1 > abs_z2) && (abs_z1 > abs_t))
  {
    xk = z[1];
    yk = S01 * e[0] + S11 * e[1] + S12 * e[2];
  }
  else if (abs_z2 > abs_t)
  {
    xk = z[2];
    yk = S02 * e[0] + S12 * e[1] + S22 * e[2];
  }
  else
  {
    xk = t;
    yk = z[0] * e[0] + z[1] * e[1] + z[2] * e[2];
  }

  // Eqn.(28)
  const double h = sqrt(xk * xk + yk * yk);
  const double h_inv = 1.0 / h;
  const double sph = xk * h_inv;
  const double cph = -yk * h_inv;

  quaternion_t q_opt;
  q_opt.q0 = cph;
  q_opt.q1 = sph * e[0];
  q_opt.q2 = sph * e[1];
  q_opt.q3 = sph * e[2];

  return q_opt;
}
