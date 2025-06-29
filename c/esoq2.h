#ifndef _ESOQ2_H_
#define _ESOQ2_H_

typedef struct
{
  double x;
  double y;
  double z;
} esoq2_vec_t;

typedef struct
{
  double q0;
  double q1;
  double q2;
  double q3;
} quaternion_t;

quaternion_t esoq2_update(const esoq2_vec_t *vb, const esoq2_vec_t *vi, double *w, const int n);

#endif // esoq2.h
