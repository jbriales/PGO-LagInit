function [V_upg,K] = metricUpgrade( V )
% [V_upg,K] = metricUpgrade( V )
% Compute matrix K such that V*K is closest to the primal feasible domain.
% 
% Note: For convenience V is treated as a column block-vector with
% (3 x k) blocks, so we store the nullspace basis V as a blkmat object.
% 
% See also blkmat

% Assert the input is a block matrix
s = size(V);
k = s(2);

if ~isa(V,'blkmat')
  n = s(1)/3;
  V = blkmat(n,1,3,k,V);
end

% Compute the quadratic form matrices
n = numel(V);
A = zeros(k^2); b = zeros(k^2,1);
for i=1:n
  W = V(i)'*V(i);
  A = A + kron(W,W);
  b = b + vec(W);
end
c = n*3;
f_vec = @(S) vec(S)'*A*vec(S) - 2*b'*vec(S) + c;

% Project quadratic form to its symmetric vectorization
% using vec(S) = P'*vecs(S)
P = build_vec2vecsProj( k );
A = P*A*P';
b = P*b;
f_vecs = @(S) vecs(S)'*A*vecs(S) - 2*b'*vecs(S) + c;

% DEBUG: Check validity of equivalent representations
% K = randn(k,k);
% sum_Frob(V,K)
% f_vec(K*K')
% f_vecs(K*K')
% keyboard

% Solve minimum for the proximity cost
s_optim = A\b;

% Undo the symmetric vectorization
S_optim = avecs(s_optim);

% Get original linear metric upgrade K so that S=K*K'
% This is done with SVD decomposition
[U,D,~] = svd(S_optim);
if k == 3
  K = U*sqrt(D);
else
  % If the rank is higher than 3, take the largest singular values
  K0 = U(:,1:3)*sqrt(D(1:3,1:3));
  % And refine the non-convex problem
  % This non-linear refinement is done with Manopt
  K = refine_fixedRank(A,b,K0);
end

% Get new linear combination
V_upg = V*K;

end

% TEMPORAL: For debug purposes
function f = sum_Frob( V,K )
% This function encodes the feasibility metric as defined in the paper

n = numel(V);
I_3 = eye(3);
f = 0;
for i=1:n
  f = f + norm( K'*V(i)'*V(i)*K - I_3, 'fro' )^2;
end

end

function [K] = refine_fixedRank(A,b,K0)
% Function with non-convex optimization of fixed-rank matrix
% for the case k>3. It uses Manopt toolbox.

k = size(K0,1);

K = symfixedrankYYfactory(k,3);

problem.M = K;
problem.cost = @cost;
problem.grad = @grad;

% checkgradient(problem);
 
% Solve.
% X0 = []; % Random for now (uninitialized)
opts = struct();
% opts.verbosity = 0;
[K, Scost, info, options] = trustregions(problem,K0,opts);


  function f = cost(Y)
    s = vecs(Y*Y');
    f = 0.5*s'*A*s - b'*s;
  end

  function g = grad(Y)
    
    % Objective part
    s = vecs(Y*Y');
    der_g_s = (A*s-b)';
    
    % Formulation part
    dervec_S_Y = ( eye(k^2) + build_vecProj(k,k) ) * kron(Y,eye(k));
    
    % Compose (chain rule)
    dervec_g_Y = der_g_s * build_vec2vecsProj(k) * dervec_S_Y;
    
    % Reshape to Y matrix
    g = reshape(dervec_g_Y, size(Y));
  end

end
