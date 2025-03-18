function [C, q] = esoq2(vb, vi, w)
%%% Second Estimator of the Optimal Quaternion (ESOQ-2)
%
% Inputs:
%   vb = Unit measurement vectors in the spacecraft body frame (3xn).
%   vi = Corresponding unit vectors known in the inertial frame (3xn).
%   w  = Non negative weight assigned to each observations (nx1).
%
% Outputs:
%   q = Optimal quaternion that transforms vectors in inertial frame to
%           vectors in body frame ([q0,q1,q2,q3]').
%   C = Optimal rotation matrix that transforms vectors in inertial
%           frame to vectors in body frame.
%
% Reference:
%   [1] Mortari - ESOQ-2, Single-Point Algorithm for Fast Optimal Spacecraft
%       Attitude Determination (1997)
%
% Rishav (2023-12-26)

  B = (vb .* repmat(w, [1,3])') * vi';
  z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];
  S = B + B'- eye(3) * trace(B); % Eqn.6
  K = [S, z; z', trace(B)];

  % Largest eigenvalue of K, lambda
  adj_BBt = det(B + B') * inv(B + B');
  adj_K = det(K) * inv(K);
  b = -2 * trace(B)^2 + trace(adj_BBt) - z' * z; % Eqn.7
  c = -trace(adj_K);
  d = det(K);
  p = (b/3)^2 + 4*d/3; % Eqn.10
  q_ = (b/3)^3 - 4*d*b/3 + c^2/2;
  u1 = 2 * sqrt(p) * cos((1/3) * acos(q_ / p^(1.5))) + b/3;
  [~, n] = size(vb);

  if (n == 2)
    g3 = sqrt(2 * sqrt(d) - b); % Eqn.15
    g4 = sqrt(-2 * sqrt(d) - b);
    lambda = (g3 + g4) / 2; % Eqn.14
  else
    g1 = sqrt(u1 - b); % Eqn.12
    g2 = -2 * sqrt(u1^2 - 4*d); % check_me
    lambda = 0.5 * (g1 + sqrt(-u1 - b - g2)); % Eqn.11
  end

  % Different possible axes
  t = trace(B) - lambda;
  S = B + B' - (trace(B) + lambda) * eye(3);
  M = t * S - z * z';
  m1 = M(:,1);
  m2 = M(:,2);
  m3 = M(:,3);
  e1 = cross(m2, m3);
  e2 = cross(m3, m1);
  e3 = cross(m1, m2);

  % Axis with greatest modulus for robustness
  if (norm(e1) > norm(e2) && norm(e1) > norm(e3))
    e = e1;
  elseif (norm(e2) > norm(e1) && norm(e2) > norm(e3))
    e = e2;
  elseif (norm(e3) > norm(e1) && norm(e3) > norm(e2))
    e = e3;
  end

  % Axis-angle for optimal quaternion
  e = e / norm(e); % Eqn.25
  x = [z; t]; % Eqn.26
  y = -[S; z'] * e; % Eqn.26
  [~, k] = max(abs(x));
  h = sqrt(x(k)^2 + y(k)^2); % Eqn.29
  sph = x(k) / h;
  cph = y(k) / h;

  q = [sph * e', cph]'; % Eqn.16
  C = quaternion_to_dcm(q);
end
