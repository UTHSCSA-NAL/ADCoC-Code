clear;
close all;

fn_txt        =   'ADCoC_Results_ARI.txt';
fd_txt        =   fopen( fullfile(fn_txt), 'at');

addpath(genpath('nmfv1_4'))

load 'sbjFea.mat'
load 'gendata\lab.mat'

minSbjFea = repmat(min(sbjFea),size(sbjFea,1),1);
maxSbjFea = repmat(max(sbjFea),size(sbjFea,1),1);
nSbjFea = (sbjFea-minSbjFea) ./ max(eps,maxSbjFea-minSbjFea);
nSbjFea(:,std(nSbjFea)==0) = [];

[numSbj, numFea] = size(nSbjFea);
NS = round( 0.8 * numSbj );

sn0 = 2*50*(1/255)^2; % A rough estimate of noise variance in the normalized data. 

%% clustering
Kcv = 4;   % number of subject sub-cluster
Krv = 8; % number of feature sub-cluster
for i = 1:length(Kcv)
    for j = 1:length(Krv)
        
        Kc = Kcv(i); Kr = Krv(j);
        fprintf(fd_txt, '\n%d Feature Clusters\n', Kr);
        
        for seedid = 7:10:37%
            
            close all;
            rand('twister',seedid);
            randid = randperm(numSbj);
            nSbjFea2 = nSbjFea(randid(1:NS), :);
            lab2 = lab(randid(1:NS));
            lab_kmeans = litekmeans(nSbjFea2, Kc,'Replicates',100);
            
            option.orthogonal = [1,1];
            option.iter = 40000; 
            option.dis = 0;
            option.residual = 1e-5;
            option.tof = 1e-5;
            
            FileName = ['SMtr_SC_' num2str(Kcv(i)) 'FC_' num2str(Krv(j)) '.mat'];
            load(FileName);
            SSMtr = sort(SMtr, 'descend');
            expec = mean(SSMtr, 2);
            Var_fea_observe = var(SSMtr,0, 2);
            
            disp('Performing Collaborative Clustering ...');
            tic;
            [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea2', Kc, Kr, option);
            toc;
            
            kLab_cc = litekmeans(Y',Kc,'Replicates',100);
            
            for inneriter = 1:2
                
                %%  Normalize Dictionaries
                sY = sqrt(sum(Y.^2,2));
                sA = sqrt(sum(A.^2));
                Y = Y ./ repmat(sY,1,size(Y,2));
                A = A ./ repmat(sA,size(A,1),1);
                S = repmat(sA',1, size(S,2)) .* S .* repmat(sY',size(S,1),1);
                
                %%  Shrink S
                [SS, inds] = sort(S(:), 'descend');
                
                scale = numSbj * numFea/(Kc*Kr);
                sn1 = sn0/inneriter^2;
                sn = sn1*scale;
                Var_fea = max(Var_fea_observe-sn, eps);
                SSr = (Var_fea.*SS + sn.*expec) ./ Var_fea_observe;
                St = S;
                St(inds) = SSr;
                Sr = reshape(St, size(S));
                
                SbjFea_new = A*Sr*Y;
                nSbjFea3 = SbjFea_new'; 
                [A,S,Y,numIter,tElapsed,finalResidual] = orthnmfrule_mod(nSbjFea3', Kc, Kr, option);
                
                if inneriter == 1
                    kLab_iter1 = litekmeans(Y',Kc,'Replicates',100);
                end
                
            end
            
            kLab = litekmeans(Y',Kc,'Replicates',100);
           
            [~,aLab] = max(A,[],2);
            
            r1 = RandIndex(lab2,kLab);
            r2 = RandIndex(lab2,lab_kmeans);
            r3 = RandIndex(lab2,kLab_cc);
            r4 = RandIndex(lab2,kLab_iter1);
                        
            %% plot clustering results
            [skLab,skInd] = sort(kLab);
            [saLab,saInd] = sort(aLab);
            figure; colormap('parula');hold on; box on;
            imagesc(nSbjFea2(skInd,saInd)); axis image; colorbar; caxis([min(nSbjFea2(:)),max(nSbjFea2(:))]); xlabel('Feature ID'); ylabel('Sbj ID');
            cpt = find(diff(skLab)~=0) + 1;
            for cpi=1:length(cpt)
                plot([0,size(nSbjFea2,2)+1],[cpt(cpi),cpt(cpi)],'r--','LineWidth',2);
            end
            
            cpa = find(diff(saLab)~=0);
            for cpi=1:length(cpa)
                plot([cpa(cpi),cpa(cpi)],[1,size(nSbjFea2,1)],'r--','LineWidth',2);
            end
            
            fprintf(fd_txt, 'KMeans:%2.3f\tCC:%2.3f\titer1:%2.3f\tADCoC:%2.3f\n', r2, r3, r4, r1);
            
        end
    end
end

disp('Finished.');
fclose all;