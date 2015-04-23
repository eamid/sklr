function K0 = gaussKernel(X,k)
% GAUSSKERNEL calculates initial Gaussian kernel matrix
% 
% K = gaussKernel(X,k)
%
% Function gaussKernel(X,k) finds the initial Gaussian kernel matrix of the
% input X. The bandwidth parameter is estimated adaptively, based on the
% k-nearest neighbor.
%
% input arguments:
% X   ----  input data matrix (N x d)
% k   ----  index of nearest neighbor for estimating sigma (default = 100)
%
% output arguments:
% K0  ----  initial kernel matrix (NxN)
%
% (C) Ehsan Amid, Aalto University

X = zscore(X); % normalize features
N = size(X,1);
if ~exist('k','var') || isempty(k)
    k = min(100,N);
end

diag_shift = 0; % apply diagonal shift on K0 for numerical stability

D = pdist2(X,X);
D_sorted = sort(D, 2, 'ascend');
sig = D_sorted(:,k); % estimate sigma (std) for each point

K0 = exp(-D.^2./(sig*sig')) + diag_shift * eye(N);
