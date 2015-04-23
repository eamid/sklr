function perms = allperms(inliers, outliers, no_reps)
% ALLPERMS returns all possible permutations
% 
% perms = allperms(inliers, outliers, reps)
%
% Function ALLPERMS returns all possible selections of two points from
% the inliers and the third point from the outliers. Number of repetitions 
% is added as the forth argument. For instance, allperms([1 2 3],[4 5], 10)
% returns:
% 
% [ 1     2     4    10
%   1     3     4    10
%   2     3     4    10
%   1     2     5    10
%   1     3     5    10
%   2     3     5    10 ]
% 
% The output of the function is used to generate inequality constraints 
% between inliers and outliers.
%
% input arguments:
% inliers      ----  set of inlier classes
% outliers     ----  set of outlier classes
% no_reps      ----  number of repetitions (default = 1)
%
% output arguments:
% perms        ----  all possible permutations
%
% (C) Ehsan Amid, Aalto University

if ~exist('no_reps','var') || isempty(no_reps)
    no_reps = 1;
end

if numel(inliers) == 1
    perms = [];
    return
else
    base = nchoosek(inliers, 2); % choose first two points from the inliers
end
last = ones(size(base,1),1) * outliers(:)'; % choose the outliers
perms = [repmat(base,length(outliers),1), last(:), no_reps * ones(numel(last),1)];
