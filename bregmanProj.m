function [K, error] = bregmanProj(K0, relative, equality, max_iter, low_rank)
% BREGMANPROJ applies Bregman projections on the kernel matrix
% 
% [K, error] = bregmanProj(K0, relative, equality, max_iter, low_rank)
%
% Function BREGMANPROJ performs iterative Bregman projections on the
% initial kernel matrix K0 until the set of inequality and equality 
% constraints are satisfied or the maximum number of iterations is reached.
%
% input arguments:
% K0         ----  initial kernel matrix (N x N)
% relative   ----  matrix of inequality constraints (N_ineq x 3)
% equality   ----  matrix of equality constraints (N_eq x 3)
% max_iter   ----  maximum number of iterations (default = 1000)
% low_rank   ----  perform low-rank approximation (default = true)
%
% output arguments:
% K          ----  transformed kernel matrix (NxN)
% error      ----  log det error D_ld(K,K0)
%
% (C) Ehsan Amid, Aalto University


if ~exist('max_iter','var') || isempty(max_iter)
    max_iter = 1000;
end

if ~exist('low_rank','var') || isempty(low_rank)
    low_rank = true;
end


no_points = size(K0,1);
if low_rank % perfrom low-rank approximation
    [Q,S,~] = svd(K0);
    eig_rat = cumsum(diag(S))/sum(diag(S));
    no_dims = find(eig_rat>=0.99,1);
    Q = Q(:,1:no_dims);
    K0 = S(1:no_dims,1:no_dims); % reduce size to (no_dims x no_dims)
    fprintf('num dims/num points = %d/%d\n',no_dims, no_points)
end
no_dims = size(K0,1);

K = K0;
[V,L,~] = svd(K0);
L = diag(L)*ones(1,no_dims);

gamma = 2;   % constant factor for inequality constraint, gamma * D(i,j)^2 <= D(i,k)^2
thr = 0.001; % threshold for equality constraints, D(i,j)^2 = D(i,k)^2 +/- thr
diag_shift = 0; % diagonal shift on K for numerical stability

unsat = false;
for iter = 1:max_iter
    no_unsat_rel = 0;
    no_unsat_eq = 0;
    rndidx = randperm(size(relative,1));
    for nr = 1:size(relative,1)
            e = full(ind2vec(relative(rndidx(nr),:),no_points));
            c1 = e * [1;-1;0]; c2 = e * [1;0;-1]; 
            d1 = e * [-1;1;0]; d2 = e * [0;1;-1]; 
            C1 = -(c2*c2') + gamma*(c1*c1'); % first constraint matrix
            C2 = -(d2*d2') + gamma*(d1*d1'); % second constraint matrix
            if low_rank % transform all matrices if low-rank approximation
                C1 = Q' * C1 * Q;
                C2 = Q' * C2 * Q;
            end
            if trace(K*C1) > 0
                no_unsat_rel = no_unsat_rel + 1;
                unsat = true;
                [U,D,~] = svd(K);
                M = U*sqrt(D);
                W = M'*C1*M;
%                 lambda = eigs(W,2); % use MATLAB eigs function
%                 a = sum(lambda)/prod(lambda)/2;
                lambda = real(eig(W));
                l1 = max(lambda);
                l2 = min(lambda);
                a = (l1+l2)/(l1*l2)/2; % Lagrange multiplier
                K = M/(eye(size(K,1))-a*W)*M' + diag_shift * eye(size(K,1)); % perform update
            end
            if trace(K*C2) > 0
                no_unsat_rel = no_unsat_rel + 1;
                unsat = true;
                [U,D,~] = svd(K);
                M = U*sqrt(D);
                W = M'*C2*M;
%                 lambda = eigs(W,2);
%                 a = sum(lambda)/prod(lambda)/2;
                lambda = real(eig(W));
                l1 = max(lambda);
                l2 = min(lambda);
                a = (l1+l2)/(l1*l2)/2;
                K = M/(eye(size(K,1))-a*W)*M' + diag_shift * eye(size(K,1));
            end
    end
    rndidx = randperm(size(equality,1));
    for ne = 1:size(equality,1)
        e = full(ind2vec(equality(rndidx(ne),:),no_points));
        for i = 1:3
            indices = setdiff(1:3,i);
            c1 = e * full(ind2vec(i,3)-ind2vec(indices(1),3));
            c2 = e * full(ind2vec(i,3)-ind2vec(indices(2),3));
            C = c1*c1' - c2*c2'; % constraint matrix
            if low_rank
                C = Q' * C * Q;
            end
            if abs(trace(K*C)) > thr
                no_unsat_eq = no_unsat_eq + 1;
                unsat = true;
                [U,D,~] = svd(K);
                M = U*sqrt(D);
                W = M'*C*M;
%                 lambda = eigs(W,2);
%                 a = sum(lambda)/prod(lambda)/2;
                lambda = real(eig(W));
                l1 = max(lambda);
                l2 = min(lambda);
                a = (l1+l2)/(l1*l2)/2;
                K = M/(eye(size(K,1))-a*W)*M' + diag_shift * eye(size(K,1));
            end
        end
    end
    fprintf('iter = %d, num of relative = %d/%d, num of equality = %d/%d\n',...
        iter, no_unsat_rel, 2 * size(relative,1), no_unsat_eq, 3 * size(equality,1));
    if ~unsat
        disp('all constraints satisfied')
        break
    else
        unsat = false;
    end
end

[U,T,~] = svd(K);
T = diag(T)*ones(1,no_dims);
error = (V'*U).^2.*(L./T' - log(L./T') - 1);
error = sum(error(:));

if low_rank 
    K = Q * K * Q'; % return back to (no_pints x no_points) matrix
end
