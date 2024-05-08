function res = logit2(beta, S, criterion, N)
  p = 2; % Dimension
  N1 = N;
  N = N1^p;

  X = cell(1, p);
  X(:) = {linspace(S(1), S(2), N1)};
  u = sortrows(combvec(X{:}).') ; 

%% 0. Initialization
tol = 10^(-4); % for finding and filtering out the points
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
  cvx_status
  res = {cvx_optval, d00, w00};
end