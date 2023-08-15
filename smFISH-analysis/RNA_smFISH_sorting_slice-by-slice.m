clear all
% For analyzing the smFISH data and calculating the amount and distribution of the RNA targets in HEK 293T and U-2 OS cells.
% Z-stack confocal imaging data of smFISH (RNA), stress granule (G3BP1-EGFP), and
% nucleus (DAPI) in 16-bit are needed
%% Import the Z-stack confocal imaging data of granule
% define the contrast of granule and smFISH (RNA) channels
cmax_granule = 2000;%
cmax_smFISH = 2000;
cmin_granule = 100
cmin_smFISH = 100;

subfolder = 'C:\Users\Luxin Peng\Desktop\134416_60x1.5 confocal-display\134416_60x1.5 confocal-display\134505_60x1.5 confocal\';
basepath = '';
pathname = [basepath subfolder];
name_granule = 'G3BP1-EGFP'; % stack (tiff file) of the granule channel should be in the subfolder, and this should be the file name
bin = 1;
mkdir([pathname 'analysis\' name_granule]);
FileTif=[pathname name_granule '.tif'];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov_granule=zeros(nImage,mImage,NumberImages,'uint16');
start_frame = 1;

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov_granule(:,:,i)=TifLink.read();
end
TifLink.close();
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
mov_granule = double(mov_granule);
mov_granule_odd = mov_granule(:,:,1:2:size(mov_granule,3));
mov_granule_even = mov_granule(:,:,2:2:size(mov_granule,3));
close(gcf)
I_odd_granule = max(mov_granule_odd,[],3);
I_even_granule = max(mov_granule_even,[],3);
I_granule = max(mov_granule,[],3);


figure()
imshow(I_odd_granule,[cmin_granule cmax_granule])
colorbar
title('Odd-projection: granule')
saveas(gcf,[pathname '\analysis\' name_granule '\max_proj._odd_granule.fig'])
saveas(gcf,[pathname '\analysis\' name_granule '\max_proj._odd_granule.png'])
close(gcf)
figure()
imshow(I_even_granule,[cmin_granule cmax_granule])
title('Even-projection: granule')
saveas(gcf,[pathname '\analysis\' name_granule '\max_proj._even_granule.fig'])
saveas(gcf,[pathname '\analysis\' name_granule '\max_proj._even_granule.png'])
close(gcf)

%% Import the Z-stack confocal imaging data of RNA (smFISH)
name_smFISH = 'APLP2'; % stack (tiff file) of the RNA channel should be in the subfolder, and this shold be the file name
% name1color = getrgb(610);
bin = 1;
basepath = '';
pathname = [basepath subfolder];
mkdir([pathname 'analysis\' name_smFISH]);
FileTif=[pathname name_smFISH '.tif'];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov_smFISH=zeros(nImage,mImage,NumberImages,'uint16');
start_frame = 1;

TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov_smFISH(:,:,i)=TifLink.read();
end
TifLink.close();
bkg = 100*power(bin,2);          % background due to camera bias (100 for bin 1x1)
mov_smFISH = double(mov_smFISH);
mov_smFISH_odd = mov_smFISH(:,:,1:2:size(mov_smFISH,3));
mov_smFISH_even = mov_smFISH(:,:,2:2:size(mov_smFISH,3));

I_odd_smFISH = max(mov_smFISH_odd,[],3);
I_even_smFISH = max(mov_smFISH_even,[],3);
figure()
imshow(I_odd_smFISH,[cmin_smFISH cmax_smFISH])
colorbar
title('Odd-projection: smFISH')
saveas(gcf,[pathname '\analysis\' name_smFISH '\max_proj._odd_smFISH.fig'])
saveas(gcf,[pathname '\analysis\' name_smFISH '\max_proj._odd_smFISH.png'])
close(gcf)
figure()
imshow(I_even_smFISH,[cmin_smFISH cmax_smFISH])
colorbar
title('Even-projection: smFISH')
saveas(gcf,[pathname '\analysis\' name_smFISH '\max_proj._even_smFISH.fig'])
saveas(gcf,[pathname '\analysis\' name_smFISH '\max_proj._even_smFISH.png'])
close(gcf)

%% Define the ROI of granule
threshold_granulemask = 2000; % a mannually-selected threshold depending on the expression
% median filter
I_granule = medfilt2(I_granule,[5,5]);
I_odd_granule = medfilt2(I_odd_granule,[5,5]);
I_even_granule = medfilt2(I_even_granule,[5,5]);
% binarization
graythreshold_granulemask = threshold_granulemask/65535;
BW_cell_granule_odd = im2bw(I_odd_granule./65535,graythreshold_granulemask);
BW_cell_granule_even = im2bw(I_even_granule./65535,graythreshold_granulemask);
BW_cell_granule = im2bw(I_granule./65535,graythreshold_granulemask);

figure()
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,3,1)
imshow(BW_cell_granule_odd)
title(subplot(1,3,1),'Granule odd: granule')
subplot(1,3,2)
imshow(BW_cell_granule_even)
title(subplot(1,3,2),'Granule even: granule')
subplot(1,3,3)
imshow(BW_cell_granule)
title(subplot(1,3,3),'Granule even: total')
saveas(gcf,[pathname '\analysis\' name_granule '\granule_bwcheck.fig'])
saveas(gcf,[pathname '\analysis\' name_granule '\granule_bwcheck.png'])

temp1 = ind2rgb(BW_cell_granule_odd*65535, pseudocolor(637));
temp2 = mat2gray(I_odd_granule,[100 cmax_granule]);
temp2 = ind2rgb(round(temp2*255), pseudocolor(500));
temp3 = ind2rgb(BW_cell_granule_even*65535, pseudocolor(637));
temp4 = mat2gray(I_even_granule,[100 cmax_granule]);
temp4 = ind2rgb(round(temp4*255), pseudocolor(500));
temp5 = ind2rgb(BW_cell_granule*65535, pseudocolor(637));
temp6 = mat2gray(I_granule,[100 cmax_granule]);
temp6 = ind2rgb(round(temp6*255), pseudocolor(500));

granule_merge_odd = temp1+temp2;
granule_merge_even = temp3+temp4;
granule_merge_total = temp5+temp6;
figure()
set(gcf,'outerposition',get(0,'screensize'));
subplot(1,3,1)
imshow(granule_merge_odd)
title(subplot(1,3,1),'Granule odd')
subplot(1,3,2)
imshow(granule_merge_even)
title(subplot(1,3,2),'Granule even')
subplot(1,3,3)
imshow(granule_merge_total)
title(subplot(1,3,3),'Granule total')
saveas(gcf,[pathname '\analysis\' name_granule '\granule_mergecheck.fig'])
saveas(gcf,[pathname '\analysis\' name_granule '\granule_mergecheck.png'])
%% select cells from the granule channel and calculation
figure()
imshow(granule_merge_odd);hold on
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
while(npts > 0)
    [xv, yv] = (getline(gca, 'closed'));
    if size(xv,1) < 3  % exit loop if only a line is drawn
        title('Selected Regions');
        hold on;
        break
    end
    currcolor = order(1+mod(colorindex,size(order,1)),:);
    plot(xv, yv, 'Linewidth',1.5,'Color',currcolor);
    text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',14);
    colorindex = colorindex+1;
    %     inpoly{nroi} = [inpolygon(x,y,xv,yv)];
    roi_points{nroi} = [xv, yv];
    nroi = nroi + 1;
end
saveas(gca,[pathname 'Selected_region.fig'])
saveas(gca,[pathname 'Selected_region.png'])
[ysize xsize] = size(BW_cell_granule_odd);
[x y] = meshgrid(1:xsize, 1:ysize);
%% Calculate the average intensity in the RNA channel of the selected cell regions slice by slice
inpoly = zeros(ysize,xsize,length(roi_points));
for n = 1:length(roi_points)
    inpoly_temp = inpolygon(x,y,roi_points{1,n}(:,1),roi_points{1,n}(:,2));
    intens(:,n) = squeeze(sum(sum(mov_smFISH.*repmat(inpoly_temp, [1, 1, size(mov_smFISH, 3)]))))/sum(inpoly_temp(:));
    inpoly(:,:,n) = inpoly_temp;
    clear inpoly_temp
end
%% Determine the smFISH localization
% The RNA-positive pixels are defined as the pixels whose intensities are
% over 2-5 times of the average intensity, where the specific times
% (threshold_smFISHmask_times) below are fixed in the analysis of
% identical RNA molecules. And the RNA-positive puncta are the connected
% components extracted from the RNA-positive pixel collections
threshold_smFISHmask_times = 3;
min_smFISH_size = 15;% The selected connected components whose areas are smaller than this value would be excluded from the RNA-positive puncta. In the papaer 12 or 15 were selected to analyzing different RNA targets
smFISH_BW = zeros(ysize,xsize,size(mov_smFISH,3),length(roi_points));
pathname_check = [pathname '\check_the_smFISH_selection\'];
mkdir(pathname_csheck)
for n = 1:size(mov_smFISH,3)
    tempslice = mov_smFISH(:,:,n);
    tempslice = medfilt2(tempslice,[5 5]);% median filter
    for cell_order = 1:length(roi_points)
        thresh_temp = 100+(threshold_smFISHmask_times*(intens(n,cell_order)-100));
        mask_temp = tempslice.*inpoly(:,:,cell_order);
        BW_temp = mask_temp > thresh_temp;
        
        [L,N] = bwlabel(BW_temp,8);  % N is the amount of the connected components
        %         L_temp = L;
        for k = 1:N
            puncta_volume(k) = length(find(L == k));
            if puncta_volume(k) <= min_smFISH_size;
                L(L == k) = 0;
            end
        end
        smFISH(:,:,n,cell_order) = L;
        % Line 211-224 are for generating the selected RNA-positive pixels in each slice 
%         figure()
%         set(gcf,'outerposition',get(0,'screensize'));
%         subplot(3,1,1)
%         imshow(mat2gray(mask_temp,[100 5*intens(n,cell_order)]))
%         title(subplot(3,1,1),['Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(3,1,2)
%         imshow(mat2gray(mask_temp,[100 5*intens(n,cell_order)])+ind2rgb(round(BW_temp*255), pseudocolor(637)))
%         title(subplot(3,1,2),['Check the smFISH selection: Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(3,1,3)
%         imshow(mat2gray(mask_temp,[100 5*intens(n,cell_order)])+ind2rgb(round(L*255), pseudocolor(637)))
%         title(subpslot(3,1,3),['Check the smFISH selection after removing small punctas: Cell ' num2str(cell_order) ', slice ' num2str(n) ', threshold = ' num2str(threshold_smFISHmask_times) 'x avg.'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_smFISHmask_times) 'x.fig'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_smFISHmask_times) 'x.png'])
%         close(gcf)
    end
end
smFISH(smFISH>0) = 1;
%% Calculation the RNA molecule amount and distribution
clear cell_size_odd granule_positive_odd smFISH_positive_odd double_positive_odd ratio_in_odd cell_size_even granule_positive_even smFISH_positive_even double_positive_even ratio_in_even
for cell_order = 1:size(roi_points,2)
    inpoly_temp = inpoly(:,:,cell_order);
    %     inpoly = roi_points{1,region};
    cell_size(cell_order) = sum(inpoly_temp(:));
    granule_positive(cell_order) = squeeze(sum(sum(BW_cell_granule.*inpoly_temp)));
    for n = 1:size(mov_smFISH,3)
    smFISH_positive(cell_order,n) = squeeze(sum(sum(smFISH(:,:,n,cell_order).*inpoly_temp)));
    double_positive(cell_order,n) = squeeze(sum(sum(BW_cell_granule.*smFISH(:,:,n,cell_order).*inpoly_temp)));
%     ratio_in_odd(cell_order,n) = double_positive_odd(cell_order)./smFISH_positive_odd(cell_order) ;
    end
end
ratio = sum(double_positive,2)./sum(smFISH_positive,2);
ratio_slice = double_positive./smFISH_positive;

figure()
for n = 1:size(intens,2)
plot(intens(:,n),'-o','Color', getrgb(400+round((n-1)/size(intens,2)*270)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
title(['Ratio in each slice']);
xlabel('Slice number')
ylabel('Average intensity of each slice')
saveas(gcf,[pathname '\analysis\Cell_intensity.fig'])
saveas(gcf,[pathname '\analysis\Cell_intensity.png'])


figure()
for n = 1:size(intens,2)
plot(ratio_slice(n,:),'-o','Color', getrgb(400+round((n-1)/size(intens,2)*270)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
title(['Ratio in each slice']);
xlabel('Slice number')
ylabel('smFISH positive ratio')
saveas(gcf,[pathname '\analysis\Ratio_distribution.fig'])
saveas(gcf,[pathname '\analysis\Ratio_distribution.png'])

figure()
for n = 1:size(intens,2)
plot(smFISH_positive(n,:),'-o','Color', getrgb(400+round((n-1)/size(intens,2)*270)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
title(['Numbers of smFISH-positive pixel in each slice']);
xlabel('Slice number')
ylabel('smFISH-positive pixel number')
saveas(gcf,[pathname '\analysis\smFISH-positive_distribution.fig'])
saveas(gcf,[pathname '\analysis\smFISH-positive_distribution.png'])

figure()
for n = 1:size(intens,2)
plot(double_positive(n,:),'-o','Color', getrgb(400+round((n-1)/size(intens,2)*270)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
title(['Numbers of double-positive pixel in each slice']);
xlabel('Slice number')
ylabel('double-positive pixel number')
saveas(gcf,[pathname '\analysis\double-positive_distribution.fig'])
saveas(gcf,[pathname '\analysis\double-positive_distribution.png'])
% 
% temptmp(:,k) = ratio_in_even;
% % k = k+1
% % end
% plot(temptmp')
%% write the results to a excel file
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],{'Region (cell) number','Cell size','Granule_positive','smFISH_positive','smFISH_in','smFISH_out','smFISH_ratio','Threshold of granule','Threshold of smFISH','Min. region size'},'Total',['A1']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],{'Region (cell) number','smFISH_positive_slice'},'smFISH_positive_slice',['A1']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],[1:1:size(roi_points,2)]','smFISH_positive_slice',['A2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],smFISH_positive,'smFISH_positive_slice',['B2']);

xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],{'Region (cell) number','double_positive_slice'},'double_positive_slice',['A1']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],[1:1:size(roi_points,2)]','double_positive_slice',['A2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],double_positive,'double_positive_slice',['B2']);

xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],{'Region (cell) number','Ratio_slice'},'Ratio_slice',['A1']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],[1:1:size(roi_points,2)]','Ratio_slice',['A2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],ratio_slice,'Ratio_slice',['B2']);


xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],[1:1:size(roi_points,2)]','Total',['A2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],cell_size','Total',['B2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],granule_positive','Total',['C2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],sum(smFISH_positive,2),'Total',['D2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],sum(double_positive,2),'Total',['E2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],sum(smFISH_positive,2)-sum(double_positive,2),'Total',['F2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],ratio,'Total',['G2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],threshold_granulemask,'Total',['H2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],threshold_smFISHmask_times,'Total',['I2']);
xlswrite([pathname 'data_analysis_threshold = ' num2str(threshold_smFISHmask_times) ', min_region = ' num2str(min_smFISH_size) '.xlsx'],min_smFISH_size,'Total',['J2']);

