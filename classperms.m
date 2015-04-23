function perms = classperms(no_classes, no_reps)
% CLASSPERMS returns all possible inter-class permutations
% 
% perms = classperms(no_clusters, no_reps)
%
% Function classperms returns all permutations of two points from one class
% and the third point from rest of the classes. The total number of classes
% is equal to no_classes. Number of repetitions is added as the forth
% argument. For instance, listperms(3,10) returns:
% 
% [ 1     1     2    10
%   1     1     3    10
%   2     2     1    10
%   2     2     3    10
%   3     3     1    10
%   3     3     2    10 ]
% 
% The output of the function is used to generate inter-class inequality constraints.
%
% input arguments:
% no_classes   ----  number of classes
% no_reps      ----  number of repetitions (default = 1)
%
% output arguments:
% perms        ----  all possible inter-class permutations
%
% (C) Ehsan Amid, Aalto University

if ~exist('no_reps','var') || isempty(no_reps)
    no_reps = 1;
end

indices = repmat(1:no_classes,no_classes,1);
perms = [indices(:) indices(:) repmat([1:no_classes]',no_classes,1)];
perms(1:no_classes+1:end,:) = []; % remove intra-class perms!

perms = [perms, no_reps * ones(size(perms,1),1)];