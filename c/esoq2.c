#include <stdio.h>
#include <math.h>

#include "esoq2.h"

static inline double det_3x3(double m[3][3])
{
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
         m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
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

  // trace(adj(B + tr(B)))
  const double tradjB =
    4 * B[1][1] * B[2][2] - (B[1][2] + B[2][1]) * (B[1][2] + B[2][1]) +
    4 * B[0][0] * B[2][2] - (B[0][2] + B[2][0]) * (B[0][2] + B[2][0]) +
    4 * B[0][0] * B[1][1] - (B[0][1] + B[1][0]) * (B[0][1] + B[1][0]);

  const double trB = B[0][0] + B[1][1] + B[2][2];
  const double trB_sq = trB * trB;
  const double z[3] = {B[1][2] - B[2][1], B[2][0] - B[0][2], B[0][1] - B[1][0]};

  // S = B + tr(B) - eye(3) * trB
  double S00 = 2 * B[0][0] - trB;
  double S11 = 2 * B[1][1] - trB;
  double S22 = 2 * B[2][2] - trB;
  const double S01 = B[0][1] + B[1][0];
  const double S02 = B[0][2] + B[2][0];
  const double S12 = B[1][2] + B[2][1];
  const double z0_sq = z[0] * z[0];
  const double z1_sq = z[1] * z[1];
  const double z2_sq = z[2] * z[2];
  const double S12_sq = S12 * S12;
  const double S01_sq = S01 * S01;
  const double S02_sq = S02 * S02;

  // trace(adj(K)) where K = [S, z; tr(z), tr(B)]
  // const double tradjK =
  //   det_3x3((double[3][3]){{S11, S12, z[1]}, {S12, S22, z[2]}, {z[1], z[2], trB}}) +
  //   det_3x3((double[3][3]){{S00, S02, z[0]}, {S02, S22, z[2]}, {z[0], z[2], trB}}) +
  //   det_3x3((double[3][3]){{S00, S01, S02}, {S01, S11, S12}, {S02, S12, S22}});

  // Determinant of 4x4 symmetric matrix K
  const double detK =
    z0_sq * S12_sq - S00 * trB * S12_sq - 2 * S01 * z[0] * z[2] * S12 +
    2 * S01 * S02 * trB * S12 + 2 * S00 * z[2] * S12 * z[1] -
    2 * S02 * z[0] * S12 * z[1] - S11 * S00 * z2_sq + S01_sq * z2_sq +
    2 * S11 * S02 * z[0] * z[2] - S11 * S02_sq * trB + S02_sq * z1_sq -
    2 * S01 * S02 * z[2] * z[1] + S11 * S00 * trB * S22 - S00 * z1_sq * S22 -
    S01_sq * trB * S22 + 2 * S01 * z[0] * z[1] * S22 - S11 * z0_sq * S22;

  const double b = -2.0f * trB_sq + tradjB - (z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
  // const double c = -tradjK;
  const double d = detK;
  // const double p = pow(b / 3.0f, 2) * 4.0f * d / 3;
  // const double q_ = pow(b / 3.0f, 3) - 4.0f * d * b / 3.0f;

  // n == 2
  const double g3 = sqrt(2 * sqrt(d) - b);
  const double g4 = sqrt(-2 * sqrt(d) - b);
  const double lambda = (g3 + g4) / 2.0f;

  const double t = trB - lambda;
  S00 -= lambda;
  S11 -= lambda;
  S22 -= lambda;

  // Eqn. 20
  const double ma = S00 * t - z0_sq;
  const double mb = S11 * t - z1_sq;
  const double mc = S22 * t - z2_sq;
  const double mx = S01 * t - z[1] * z[0];
  const double my = S02 * t - z[0] * z[2];
  const double mz = S12 * t - z[1] * z[2];

  double e[3] = {mb * mc - mz * mz, ma * mc - my * my, ma * mb - mx * mx};
  double e_abs[3] = {fabs(e[0]), fabs(e[1]), fabs(e[2])};

  // Choose optimal principal axis
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

  const double norm_e = sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  e[0] = e[0] / norm_e;
  e[1] = e[1] / norm_e;
  e[2] = e[2] / norm_e;

  const float abs_z0 = fabs(z[0]);
  const float abs_z1 = fabs(z[1]);
  const float abs_z2 = fabs(z[2]);
  const float abs_t = fabs(t);

  double xk = 0.0;
  double yk = 0.0;

  if((abs_z0 > abs_z1) && (abs_z0 > abs_z2) && (abs_z0 > abs_t))
	{
    xk = z[0];
		yk = S00 * e[0] + S01 * e[1] + S02 * e[2];
	}
  else if((abs_z1 > abs_z2) && (abs_z1 > abs_t))
	{
    xk = z[1];
		yk = S01 * e[0] + S11 * e[1] + S12 * e[2];
	}
  else if(abs_z2 > abs_t)
	{
    xk = z[2];
		yk = S02 * e[0] + S12 * e[1] + S22 * e[2];
	}
  else
	{
    xk = t;
		yk = z[0] * e[0] + z[1] * e[1] + z[2] * e[2];
	}

  const double  h = sqrt(xk * xk + yk * yk);
  const double sph = xk / h;
  const double cph = -yk / h;

  quaternion_t q;
  q.q0 = cph;
  q.q1 = sph * e[0];
  q.q2 = sph * e[1];
  q.q3 = sph * e[2];

  return q;
}
