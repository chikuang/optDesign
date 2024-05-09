%% Author: Chi-Kuang Yeh
% Sep 11 / 2017
% C optimality criterion
% Optimal designs for regression models using 2nd order LSE

% [a,b,c] = c_opt_SLSE(101, 0, [0.05, 0.5]', [0,180]',@peleg, [0,1])
%% function itself
function [del, ANS, error] = c_opt_SLSE(N, t, theta, range, fun, c)
  %% initialization
  u = linspace(range(1), range(2), N); %discretized equally spaced space
  w = zeros(N, 1); del = 0; n = length(theta);
  g1 = zeros(n, 1); G2 = zeros(n); C = [0, c]; one_vec = ones(N, 1); zero_vec = zeros(N, 1);
  %C is only used in C-optimality, for paramter combination

  %% cvx part
  cvx_begin
  cvx_precision best
  variables w(N, 1) del(1) %design variable and upper bound
  minimize del(1)
  subject to
  % constructing the B matrix
  for i = 1:N
    f = fun(u(i), theta);
    g1 = g1 + w(i) * f;
    G2 = G2 + w(i) * f * f';
  end

  B = [1, sqrt(t) * g1'; sqrt(t) * g1, G2];
  % the three constrains
  matrix_frac(C', B) <= del;
  %-w <= zeros(length(w),1);
  -w <= zero_vec;
  one_vec' * w == 1;
  cvx_end

  %% checking condition
  BI = inv(B);
  phi_C = zeros(N, 1);

  for i = 1:N
    f = fun(u(i), theta);
    I = [1, sqrt(t) * f'; sqrt(t) * f, f * f'];
    phi_C(i) = trace(I * BI * C' * C * BI);
  end

  % update the error
  term = C * BI * C';
  error = max(phi_C - term);
  %% out the result
  kk = find(w > 1e-2); % because the computer does not have exact zero
  ANS = [u(kk); w(kk)']; % return the answer

  %% plot
  % first, we increase the graphing domain
  new_range = [0; 0]; add_dist = (range(2) - range(1)) / 20;
  new_range(1) = range(1) - add_dist;
  new_range(2) = range(2) + add_dist;

  % first plot, plot out the design points
  figure
  stem(u, w, 'kd')
  xlim(new_range); % increase the domain, so the points wont be on the edge
  xlabel('design space', 'FontSize', 16) % x-axis label
  ylabel('weight', 'FontSize', 16) % y-axis label
  title('Discretized weight distribution', 'FontSize', 20)

  %directional derivative
  fx = @(x) fun(x, theta); sqt = sqrt(t);
  ff = @(x) trace([1 sqrt(t) * fx(x)'; sqrt(t) * fx(x) fx(x) * fx(x)'] * BI * C' * C * ...
    BI) - term;

  figure
  h1 = plot(u, phi_C - term * ones(N, 1), '+'); %discretized
  hold on
  h2 = fplot(ff, [range(1) range(2)], '-'); %function
  hold on
  line([range(1), range(2)], [0, 0], 'Color', 'red');
  hold on
  h3 = plot(u(kk), zeros(1, length(kk)), 'pg'); %supporting points
  legend([h1 h2 h3], 'Discretized', 'd(x,\theta)', 'Supporting point');
  xlabel('design space', 'FontSize', 16); % x-axis label
  ylabel('Directional Derivative', 'FontSize', 16); % y-axis label
  hold off
end