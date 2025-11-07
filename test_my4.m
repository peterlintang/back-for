function [Ainv, b, r] = test_my4(x, y, z)
  n = length(x);
  X = [2*x, 2*y, 2*z, ones(n, 1)];
  Y = x.^2 + y.^2 + z.^2;
  XTY = X' * Y;
  XTX = X' * X;
  v = -XTX \ XTY;

  u = zeros(10, 1);
  u(1:3) = ones(size(u(1:3)));
  u(7:9) = v(1:3);
  u(10) = v(4);

  [A, Ainv, b] = get_params(u);
  r = cal_r(u, A, b);
end

function [A, Ainv, b] = get_params(v)
  A = [v(1), v(6), v(5); v(6), v(2), v(4); v(5), v(4), v(3)];
  k = [v(7); v(8); v(9)];

  Ainv = sqrtm(A ./ nthroot(det(A), 3));

  b = -A \ k;
end
function [r] = cal_r(v, A, b)
  terms =  A(1,1) * b(1)^2 + 2 * A(1,2) * b(1) * b(2) ...
         + A(2,2) * b(2)^2 + 2 * A(1,3) * b(1) * b(3) ...
         + A(3,3) * b(3)^2 + 2 * A(2,3) * b(2) * b(3) - v(10);
  r = sqrt(abs(terms)) / sqrt(nthroot(det(A), 3));
end
