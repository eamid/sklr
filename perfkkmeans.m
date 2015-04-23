function [IDX,err_min] = perfkkmeans(K, no_clusters, no_reps)
% PERFKKMEANS repeats kernel k-means algorithm
%
% [IDX,err_min] = perfkkmeans(K, no_clusters, no_reps)
%
% Function PERFKKMEANS repeats kernel k-means algorithm no_reps times. 
% 
% input arguments:
% K            ----  kernel matrix (N x N)
% no_clusters  ----  number of clusters
% no_reps      ----  number of repetitions
%
% output arguments:
% IDX          ----  cluster indices (N x 1)
% error_min    ----  minimum error over all repetitions
%
% (C) Ehsan Amid, Aalto University

err = zeros(no_reps,1);
idx_temp = zeros(size(K,1),no_reps);

for tt = 1:no_reps
    fprintf('iteration %3d out of %3d ... ',tt,no_reps)
    [idx_temp(:,tt), err(tt)] = kernelkmeans(K, no_clusters, 1, 0);
    fprintf('done\n')
end

[err_min,idx_min] = min(err);
IDX = idx_temp(:,idx_min);
