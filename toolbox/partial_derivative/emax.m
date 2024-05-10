% Author : Chi-Kuang Yeh
% Date: 2017 Sep 20th

% Purpose: To calculate the gradient of the function in EMAX model
% Theta has dimension 3x1

function ANS = emax(xi, theta)
  S = 1 / (theta(2) + xi^theta(3));
  X_t3 = 1 / S - theta(2);
  X = (1 / S - theta(2))^(1 / theta(3));
  par_1 = X_t3 * S;
  par_2 = -theta(1) * X * S^2;
  par_3 = theta(1) * theta(2) * X / theta(3) * log(X) * S^2;
  ANS = [par_1; par_2; par_3];
end