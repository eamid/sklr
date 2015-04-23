function [IDX, error, centers] = kernelkmeans(K, no_clusters, init, disp_output)
% KERNELKMEANS applies kernel k-means algorithm
%
% [IDX, error] = kernelkmeans(K, no_clusters, init, disp_output)
%
% Function KERNELKMEANS applies the kernel k-means algorithm on kernel
% matrix K. The output indices are stores in IDX.
% 
% input arguments:
% K            ----  kernel matrix (N x N)
% no_clusters  ----  number of clusters
% init         ----  farthest neighbor search for initialization (default = true)
% disp_output  ----  display error in each iteration (default = true)
%
% output arguments:
% IDX          ----  cluster indices (N x 1)
% error        ----  sum of squared errors
%
% (C) Ehsan Amid, Aalto University

N = size(K,1);
IDX = zeros(N,1);
thr = 1e-6; % threshold for stopping

if ~exist('disp_output', 'var') || isempty(disp_output)
    disp_output = true;
end

% initialize

if ~exist('init', 'var') || isempty(init)
    init = true;
end

if init
    if disp_output
        disp('initalize...')
    end
    centers = zeros(1,no_clusters);
    centers(1) = randi(N);
    IDX(centers(1)) = 1;
    if no_clusters > 1
        D = bsxfun(@plus, diag(K), diag(K)') - 2 * K;
        D(D < 0) = 0;
        D(1:N+1:end) = 0;
        [~,sortlist] = sort(D,2,'ascend');
        for cu = 2:no_clusters
            for last = cu-1:N-1
                found = true;
                isfar = true(1,last);
                for ck = 2:cu-1
                    hits = ismember(sortlist(centers(1),end-last+1:end), sortlist(centers(ck),end-last+1:end));
                    if ~any(hits)
                        found = false;
                        break
                    else
                        isfar = isfar & hits;
                    end
                end
                if found && any(isfar)
                    break
                end
            end
            sublist = sortlist(centers(1),end-last+1:end);
            elem = find(isfar);
            centers(cu) = sublist(elem(end));
            IDX(centers(cu)) = cu;
        end
    end
    if disp_output
        disp('end initialize.')
    end
else
    centers = randperm(N);
    centers = centers(1:no_clusters);
    IDX(centers) = 1:no_clusters;
end

old_error = Inf;
dist = zeros(N,no_clusters);
sum_diag = sum(diag(K));

if disp_output
    disp('start clustering')
end
iter = 1;
while 1
    for cu = 1:no_clusters
        hit = IDX == cu;
        list = bsxfun(@times, hit, hit');
        const = sum(K(list == 1))/sum(hit)^2;
        dist(:,cu) = -2 * sum(K(:,hit),2)/sum(hit) + const;
    end
    [error, IDX] = min(dist,[],2);
    error = sum(error) + sum_diag;
    if abs(old_error - error)/old_error < thr
        break
    else
        old_error = error;
    end
    if disp_output
        fprintf('iteration %d, error = %f\n', iter, log(error));
    end
    iter = iter + 1;
end
if disp_output
    disp('end clustering')
end