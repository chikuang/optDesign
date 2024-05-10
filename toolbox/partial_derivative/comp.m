% Author : Chi-Kuang Yeh 
% Date: 2017 July 10th

% Purpose: To calculate the gradient of the function in 
%   Compartmental model with n/2 compartments
% Theta has dimension nx1

% [d, a, e] = D_opt(101, 0.3, [1, 0.1,1, 0.6, 1, 2.3, 1, 5.5]', [0;1], @comp)
function ANS = comp(xi, theta)
  ANS = zeros(length(theta), 1);
  n = length(theta);
  count  = 0 ;
  i = 1;
  while( count < (n/2) )
    ANS(i, 1) = exp(-theta(i+1) * xi);
    ANS(i+1, 1) =  -theta(i) * xi * exp(-theta(i+1) * xi);
    i = i + 2;
    count = count + 1;
  end
end