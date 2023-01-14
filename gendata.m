bs = 37;
ori = 20*ones(bs, bs);
r1 = 2;
r2 = 3;
r3 = 7;

r4 = 5;
r5 = 11;

t1 = 50;
t2 = 200;
t3 = 250;
vnoise = sqrt(50);
lab = [];
sbjFea = [];
sbjFea_clean = [];
stride = 4;
stride2 = 3;

k = 0;
for c1 = r4+3:stride:bs-r4-3
    for c2 = r4+3:stride:bs-r4-3
        k = k+1;
        ori_tem = ori;
        ori_tem(c1-r4:c1+r4, c2-r4:c2+r4) = t2;
        lab = [lab;1];
        sbjFea_clean = [sbjFea_clean; ori_tem(:)'];
        ori_tem = ori_tem + vnoise * randn(size(ori_tem));
        imfea = imresize(ori_tem,0.5);
        sbjFea = [sbjFea; imfea(:)'];
        imagesc(ori_tem./255);axis off; caxis([0,1])
        saveas(gcf,['gendata\' num2str(r4) '_' num2str(k) '.png']);
    end
end

k = 0;
for c1 = r3+stride2:stride2:bs-r3-stride2
    for c2 = r3+stride2:stride2:bs-r3-stride2
        k = k+1;
        ori_tem = ori;
        [xx,yy] = ndgrid((1:bs)-c1,(1:bs)-c2);
        mask = uint8((xx.^2 + yy.^2)<r3^2);
        ori_tem(mask==1) = t2;
        lab = [lab;2];
        sbjFea_clean = [sbjFea_clean; ori_tem(:)'];
        ori_tem = ori_tem + vnoise * randn(size(ori_tem));
        imfea = imresize(ori_tem,0.5);
        sbjFea = [sbjFea; imfea(:)'];
        
        imagesc(ori_tem./255);axis off; caxis([0,1])
        saveas(gcf,['gendata\' num2str(r3) '_' num2str(k) '.png']);
    end
end


%%
k = 0;
for c1 = r2+stride:stride:bs-r2-stride
    for c2 = r2+stride:stride:bs-r2-stride
        k = k+1;
        ori_tem = ori;
        ori_tem(c1-r2:c1+r2, c2-r2:c2+r2) = t2;
        lab = [lab;3];
        sbjFea_clean = [sbjFea_clean; ori_tem(:)'];
        ori_tem = ori_tem + vnoise * randn(size(ori_tem));
        imfea = imresize(ori_tem,0.5);
        sbjFea = [sbjFea; imfea(:)'];
        imagesc(ori_tem./255);axis off; caxis([0,1])
        saveas(gcf,['gendata\' num2str(r2) '_' num2str(k) '.png']);
    end
end

k = 0;
stride2 = 2;
for c1 = r5+stride2:stride2:bs-r5-stride2
    for c2 = r5+stride2:stride2:bs-r5-stride2
        k = k+1;
        ori_tem = ori;
        [xx,yy] = ndgrid((1:bs)-c1,(1:bs)-c2);
        mask = uint8((xx.^2 + yy.^2)<r5^2);
        ori_tem(mask==1) = t2;
        lab = [lab;4];
        sbjFea_clean = [sbjFea_clean; ori_tem(:)'];
        ori_tem = ori_tem + vnoise * randn(size(ori_tem));
        imfea = imresize(ori_tem,0.5);
        sbjFea = [sbjFea; imfea(:)'];
        
        imagesc(ori_tem./255);axis off; caxis([0,1])
        saveas(gcf,['gendata\' num2str(r5) '_' num2str(k) '.png']);
    end
end
%%

save('gendata\lab.mat','lab');
save('sbjFea.mat','sbjFea');
save('sbjFea_clean.mat','sbjFea_clean');