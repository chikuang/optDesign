%% Compute D-optimal designs for logistic regression with interaction

%% Reference paper: 
% Haines, L. M., & Kabera, G. M. (2018). D-optimal designs 
% for the two-variable binary logistic regression model with interaction. Journal of Statistical Planning and Inference, 193, 136-150.
 

criterion = "D";
% beta = [1, 2, 2, 0.2]';
beta = [-3, 4, 6, 1]' ; %Example 4.2 (b)
% beta = [-2.2054, 13.5803, 2.2547, 1.6262]' ; %A real example
S1 = [0, 1]; 
S2 = [0, 1]; 
p = 2; % Dimension
Nsim = 100;
N1 = 51; % Number of design p       oints of each dimension
N = N1^p;
n = [15]';

%% Generate multi-dimensional space (using some Matlab tricks)
X = cell(1, p);
X(:) = {linspace(S1(1), S1(2), N1)};
u = sortrows(combvec(X{:}).') ; 

%% 0. Initialization
tol = 10^(-4); % for finding and filtering out the points 
tol_annealing = 1E-40;
q = length(beta);

%The following vectors and matrices are used in the information
%matrices below.

%% Step 2, calculate ethe optimal design of in S_t
cvx_begin quiet
cvx_precision best
  % variables wk(Np, 1) del(1) 
  variable w(N)
  expression M(q,q); 
  % minimize del(1)
  % subject to
  for i = 1:N
      xx = u(i,:);
      rx = [1, xx, xx(1) * xx(2)]';
      Gamma = exp(beta' * rx)/(1+exp(beta' * rx))^2;
      M = M + w(i) * Gamma * (rx * rx');
  end
  
   if criterion == "D"
    minimize (-log(det_rootn(M)));
  elseif criterion == "A"
    minimize( trace_inv(M) );   %A-opt
   elseif criterion == "E"
     minimize(-lambda_min(M))
  else
    fprintf('Does not run.');
  end
  ones(N, 1)' * w == 1;
  -w <= zeros(N, 1);
cvx_end
        

% organizing the result
kk = find(w > tol); % because the computer does not have exact zero
design_app = [u(kk,:), w(kk)]; %optimal design 
d00 = design_app(:, 1:end-1); % support points
w00 = design_app(:, end); % optimal weight
L00 = cvx_optval; % optimal objective value
cvx_status