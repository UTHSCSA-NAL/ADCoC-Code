% clc;
clear;
close all;

addpath(genpath('nmfv1_4'))

load 'sbjFea.mat'
% sbjFea(sbjFea<0) = 0;
minSbjFea = repmat(min(sbjFea),size(sbjFea,1),1);
maxSbjFea = repmat(max(sbjFea),size(sbjFea,1),1);
nSbjFea = (sbjFea-minSbjFea) ./ max(eps,maxSbjFea-minSbjFea);
nSbjFea(:,std(nSbjFea)==0) = [];
% nSbjFea = sbjFea./255;

[numSbj, numFea] = size(nSbjFea);
NS = round( 0.8 * numSbj ) ;

%% clustering
Kcv = 4;   % number of subject sub-cluster
Krv = 8; % number of feature sub-cluster

for i = 1:length(Kcv)
    for j = 1:length(Krv)
        
        SMtr = [];
        Kc = Kcv(i); Kr = Krv(j);
        
        for k = 1:100
            
            randid = randperm(numSbj);
            nSbjFea2 = nSbjFea(randid(1:NS), :);
            
            close all;
                        
            lambda = 2^-9.5;
            option.orthogonal = [1,1];
            option.iter = 300000;%40000
            option.dis = 0;
            option.residual = 1e-5;
            option.tof = 1e-5;
            
            for inneriter = 1:1
                
                disp('Performing Collaborative Clustering ...');
                tic;
                %             rand('twister',7);
                [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea2', Kc, Kr, option);
                toc;
                
                %%  Normalize Dictionaries
                sY = sqrt(sum(Y.^2,2));
                sA = sqrt(sum(A.^2));
                Y = Y ./ repmat(sY,1,size(Y,2));
                A = A ./ repmat(sA,size(A,1),1);
                S = repmat(sA',1, size(S,2)) .* S .* repmat(sY',size(S,1),1);
                
                %%  
                [SS, inds] = sort(S(:), 'descend');
                SMtr = [SMtr SS];
                
            end

        end
        
        save(['SMtr_SC_' num2str(Kc) '_FC_' num2str(Kr) '.mat'],'SMtr');

    end
end

disp('Finished.');
