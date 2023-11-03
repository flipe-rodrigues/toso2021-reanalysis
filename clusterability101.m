close all;
clear;
clc;
%%
N = 1e2;
mu1 = randn(2,1);
mu2 = randn(2,1)*2;
mu3 = randn(2,1)*4;
s1 = rand(2);
s2 = rand(2);
s3 = rand(2);
sigma1 = s1 * s1';
sigma2 = s2 * s2';
sigma3 = s3 * s3';
sigma1 = eye(2);
sigma2 = eye(2);
sigma3 = eye(2);
R1 = mvnrnd(mu1,sigma1,N);
R2 = mvnrnd(mu2,sigma2,N);
R3 = mvnrnd(mu3,sigma3,N);
R = [R1;R2;R3];
L = [zeros(N,1);ones(N,1);ones(N,1)*2];
figure;
gscatter(R(:,1),R(:,2),L);
%%
D = pdist(R,'euclidean');
figure;
histogram(D,40);
%%
[W,S] = pca(R);
figure;
histogram(S(:,1),40);
% figure;
% gscatter(S(:,1),S(:,2),L);
%% silhouette score
k_min = 2;
k_max = 6;
k_idcs = k_min : k_max;
k_n = numel(k_idcs);
% preallocation
silhouette_score = nan(N,k_n);
distortion = nan(N,k_n);
for kk = 1 : k_n
    [cluster_ids,C,sumd,D] = kmeans(R,k_idcs(kk));
    %     gm = fitgmdist(eval_mat,k_test);
    %     ids = cluster(gm,R);
    figure;
    gscatter(R(:,1),R(:,2),cluster_ids);
    for nn = 1 : N
        intra_cluster_flags = cluster_ids == cluster_ids(nn);
        dists2centroids = pdist2(R(nn,:),C);
        dists2centroids(cluster_ids(nn)) = nan;
        [~,nearest_id] = min(dists2centroids);
        nearest_neighbour_flags = cluster_ids == nearest_id;
        d = mean(pdist2(R(nn,:),R(intra_cluster_flags,:)));
        dp = mean(pdist2(R(nn,:),R(nearest_neighbour_flags,:)));
        silhouette_score(nn,kk) = (dp - d) / max(dp,d);
        distortion(nn,kk) = sumd(cluster_ids(nn));
    end
end
figure;
hold on;
plot(k_min:k_max,silhouette_score,'.k')
plot(k_min:k_max,mean(silhouette_score,1),...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor','w',...
    'linewidth',1.5)
figure;
hold on;
plot(k_min:k_max,distortion,'.k')
plot(k_min:k_max,mean(distortion,1),...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor','w',...
    'linewidth',1.5)
%%
clusterability = evalclusters(R,'kmeans','gap',...
    'referencedistribution','uniform',...
    'klist',1:k_max);
figure;
hold on;
plot(clusterability)
plot(clusterability.OptimalK,...
    clusterability.CriterionValues(clusterability.OptimalK),...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor','w',...
    'linewidth',1.5)