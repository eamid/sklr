function [relative, equality] = genConstraint(labels, rel_class, eql_class)
% GENCONSTRAINT generates random equality and inequality constraints
% 
% [relative, equality] = genConstraint(labels, rel_class, eql_class)
%
% Function GENCONSTRAINT generates random equality and inequality
% constraints based on the EQL_CLASS and REL_CLASS labels, respectively.
% Variable LABELS is the labels of the datapoints. EQL_CLASS is a (N_eq x 4)
% matrix where each row is of the form [i,j,k,no_reps], in which i, j, and
% k are the labels of the equal classes. REL_CLASS is of the form
% [i,j,k,no_reps], in which k is the outlier class among i, j, and k. Use
% functions classperms, allperms and eqlperms to generate class indices.
% 
%
% input arguments:
% labels       ----  labels of the datapoints
% rel_class    ----  labels of the classes in inequality constraints
% eql_class    ----  labels of the classes in equality constraints
%
% output arguments:
% relative     ----  inequality constraints (N_ineq x 3), N_ineq = sum(rel_class(:,3))
% equality     ----  equality constraints (N_eq x 3), N_eq = sum(eql_class(:,3))
%
% (C) Ehsan Amid, Aalto University

if ~exist('rel_class', 'var') || isempty(rel_class)
    rel_class = [];
    relative = [];
else
    relative = zeros(sum(rel_class(:,4)),3);
end

if ~exist('eql_class', 'var') || isempty(eql_class)
    eql_class = [];
    equality = [];
else
    equality = zeros(sum(eql_class(:,4)),3);
end

cnt = 0;
for i = 1:size(rel_class,1)
    triplet = rel_class(i,1:3);
    n = rel_class(i,4);
    id1 = find(labels == triplet(1)); l1 = length(id1);
    id2 = find(labels == triplet(2)); l2 = length(id2);
    id3 = find(labels == triplet(3)); l3 = length(id3);
    const = [id1(randi(l1,n,1)), id2(randi(l2,n,1)), id3(randi(l3,n,1))];
    relative((1:n)+cnt,:) = const;
    cnt = cnt + n;
end

cnt = 0;
for i = 1:size(eql_class,1)
    triplet = eql_class(i,1:3);
    n = eql_class(i,4);
    id1 = find(labels == triplet(1)); l1 = length(id1);
    id2 = find(labels == triplet(2)); l2 = length(id2);
    id3 = find(labels == triplet(3)); l3 = length(id3);
    const = [id1(randi(l1,n,1)), id2(randi(l2,n,1)), id3(randi(l3,n,1))];
    equality((1:n)+cnt,:) = const;
    cnt = cnt + n;
end
