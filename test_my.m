function [Ainv, b] = test_my(x, y, z)
  % Design matrix
  D = [x.^2, y.^2, z.^2, 2*y.*z, 2*x.*z, 2*x.*y, 2*x, 2*y, 2*z, ones(size(x))]';

  % Inverse of constraint matrix with k = 4, Eqn(7)
  invC = [0,   0.5, 0.5, 0,    0,    0;
          0.5, 0,   0.5, 0,    0,    0;
          0.5, 0.5, 0,   0,    0,    0;
          0,   0,   0,  -0.25, 0,    0;
          0,   0,   0,   0,   -0.25, 0;
          0,   0,   0,   0,    0,   -0.25];

  % Eqn(11)
  S = D * D';
  S11 = S(1:6, 1:6);   % 6X6
  S12 = S(1:6, 7:10);  % 6X4
  S22 = S(7:10, 7:10); % 4X4
  Sx = pinv(S22) * S12';

  % Eqn(14) and Eqn(15)
  M = invC * (S11 - S12 * Sx);
  [V, L] = eig(M);

  % Index of the 'only' positive eigenvalue.
  [max_ridx, ~] = max(L);
  [~, max_cidx] = max(max_ridx);

  % Eigenvector corresponding to max_cidx
  u1 = V(:, max_cidx);
  u2 = -Sx * u1;
  u = [u1', u2']';

  % Ellipsoid in matrix form: Ax + k = 0
  det_A = det([u(1), u(6), u(5); u(6), u(2), u(4); u(5), u(4), u(3)]);

  % Ensure positive definiteness
  if det_A < 0
    u = -u;
  end
  [Ainv, b] = get_calibration_params(u);
  Ainv = real(Ainv);
end

% Ellipsoid coefficients to calibration params
function [Ainv, b] = get_calibration_params(v)
  % Ellipsoid in matrix form: Ax + k = 0
  A = [v(1), v(6), v(5); v(6), v(2), v(4); v(5), v(4), v(3)];
  k = [v(7); v(8); v(9)];

  % Soft-iron correction matrix
  Ainv = sqrtm(A ./ nthroot(det(A), 3));

  % Hard-iron bias vector
  b = -A \ k;
end
