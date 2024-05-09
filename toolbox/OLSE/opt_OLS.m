% [a, b] = opt_OLS("D", 1001, [0, 4], [1,1]', @mm)
function [optval, design] = opt_OLS(criterion, N, S, beta, fun)
a = S(1); b = S(2);
f = fun;
tol = 1E-4; % for finding and filtering out the points
q = length(beta);
% q = p+1; % how many beta's (degree + 1 intercept term)
u = linspace(a, b, N); %equally spaced N points in [a,b]

cvx_begin
cvx_precision best
variable w(1, N);
expression A(q, q);

% here we compute the information matrix at
for j = 1:N
  f1 = f(u(j), beta);
  A = A + w(j) * (f1 * f1');
end

if criterion == "D"
  minimize(-log(det_rootn(A)));
elseif criterion == "A"
  minimize(trace_inv(A)); %A-opt
else
  fprintf('Does not run.');
end
0 <= w <= 1;
sum(w) == 1;
cvx_end

% organizing the result
design = [u(find(w(1, :) > tol)); w(find(w(1, :) > tol))]; %optimal design
optval = cvx_optval; % optimal objective value
end
