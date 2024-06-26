% [opt_val, opt_design, error] = opt_SLSE("D", 101, 0, [0.05, 0.5]', [0,180]',@peleg)
%% Author: Chi-Kuang Yeh
% May 19 / 2017
% A optimality criterion
% Optimal designs for regression models using 2nd order LSE

%% function itself
function [del, ANS, error] = opt_SLSE(criterion, N, t, theta, range, fun)
  %% initialization
  u = linspace(range(1), range(2), N); %discretized equally spaced space
  w = zeros(N, 1); n = length(theta); del = 0;

  if criterion == 'A'
    C = [zeros(n, 1) eye(n)]'; % each column is a ci which we used later
  elseif criterion == "c"
    C =  [0, c];
  end
  
  g1 = zeros(n, 1); G2 = zeros(n); obj_val = 0; one_vec = ones(N, 1); zero_vec = zeros(N, 1);

  %% cvx part
  cvx_begin quiet
  cvx_precision best
  variables w(N, 1) del(1) %design variable
  minimize del(1)
  subject to

  for i = 1:N
    f = fun(u(i), theta);
    g1 = g1 + w(i) * f;
    G2 = G2 + w(i) * f * f';
  end

  B = [1, sqrt(t) * g1'; sqrt(t) * g1, G2];
  % three constains

 if criterion == "A"
  for k = 1:n
    obj_val = obj_val + matrix_frac(C(:, k), B);
  end
 elseif criterion == "c"
   obj_val = matrix_frac(C', B);
 elseif criterion == "D"
   obj_val = -log_det(B);
 end

  % the three constrains
  obj_val <= del;
  %-w <= zeros(length(w),1);
  -w <= zero_vec;
  one_vec' * w == 1;
  cvx_end

  %% manage the outputs
  kk = find(w > 1e-4); % because the computer does not have exact zero
  ANS = [u(kk); w(kk)'];

  %% checking condition, from M-A
  % prepare the variables
  BI = inv(B);
  phi = zeros(N, 1);

  if criterion == "A"
    C = blkdiag(0, eye(n));
  elseif criterion == "c"
    C = [0, c];
  end
  for i = 1:N
    f = fun(u(i), theta);
    I = [1 sqrt(t) * f'; sqrt(t) * f f * f'];
    if criterion == "A" | criterion == "c"
      phi(i) = trace(I * BI * C' * C * BI);
    elseif criterion == "D"
      phi(i) = trace(B \ I);
    end
  end

  % update the error
  if criterion == "A"
    term = trace(C * BI * C');
    error = max(phi - term);
  elseif criterion == "c"
     term = C * BI * C';
     error = max(phi - term);
  elseif criterion == "D"
    q = (n + 1) * ones(N, 1);
    term = q;
    error = max(phi - q);
  end
  
  %% plots
  % first, we increase the graphing domain
  new_range = [0; 0]; add_dist = (range(2) - range(1)) / 20;
  new_range(1) = range(1) - add_dist;
  new_range(2) = range(2) + add_dist;

  % first plot, plot out the design points
  figure % new figure window
  stem(u, w, 'kd');
  xlim(new_range); % increase the domain, so the points wont be on the edge
  ylim([0, 1]);
  xlabel('design space', 'FontSize', 16) % x-axis label
  ylabel('weight', 'FontSize', 16) % y-axis label
  title('Discretized weight distribution', 'FontSize', 20)

  %directional derivative plot
   fx = @(x) fun(x, theta);
  if criterion == "A"
    ff = @(x) trace([1 sqrt(t) * fx(x)'; sqrt(t) * fx(x) fx(x) * fx(x)'] * BI * C' * C * BI) - term;
  elseif criterion == "c"
      ff = @(x) trace([1 sqrt(t) * fx(x)'; sqrt(t) * fx(x) fx(x) * fx(x)'] * BI * C' * C * ...
    BI) - term;
  elseif criterion == "D"
    ff = @(x) trace(B \ [1 sqrt(t) * fx(x)'; sqrt(t)* fx(x) fx(x) * fx(x)']) - (n + 1);
  end
  mini = min(phi - term * 1);
  figure
  h1 = plot(u, phi - term * 1, '+'); %discretized
  xlim(new_range);
  ylim([mini + mini / 10, 1]);
  hold on
  y = zeros(size(u));

  for i = 1:length(u)
    y(i) = ff(u(i));
  end

  h2 = plot(u, y, '-'); %function
  %h2 = fplot(ff,range','-'); %not taking vector, not good
  hold on
  line(new_range, [0, 0], 'Color', 'blue', 'LineStyle', '--');
  hold on
  h3 = plot(u(kk), zeros(1, length(kk)), 'pg'); %supporting points
  legend([h1 h2 h3], 'Discretized', 'd(x,\theta)', 'Supporting point');
  xlabel('design space', 'FontSize', 16); % x-axis label
  ylabel('Directional Derivative', 'FontSize', 16); % y-axis label
  hold off
end