function [Ainv, b] = test_my7(x, y, z)
  n = length(x);
  D = [x.^2, y.^2, z.^2, 2*x, 2*y, 2*z, ones(n, 1)];
  S = D' * D;

  % Eigenvector to smallest eigenvalue
  [V, L] = eig(S);
  [~, idx] = min(diag(L));
  v = V(:, idx);

  % Extract ellipsoid parameters
  u = zeros(10, 1);
  u(1:3) = v(1:3); % a, b, c
  u(7:9) = v(4:6); % p, q, r
  u(10) = v(7); % d

  % Ensure positive definiteness of the ellipsoid matrix
  det_A = det([u(1), u(6), u(5); u(6), u(2), u(4); u(5), u(4), u(3)]);

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
