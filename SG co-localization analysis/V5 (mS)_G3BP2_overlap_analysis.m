clear all
%% Loading the data from original tiff stacks
% The filename of three stacks should be input at line 16, 23, and 30 noted
% as "Keyword 1", "Keyword 2", and "Keyword 3", which are pointed to DAPI, V5
% (mS) and G3BP2-EGFP tiff stack in this script, respectively.
% These function scripts are needed: autorange.m; clicky_v2.m; getrgb.m; pseudocolor.m;
% RScanDir.m; 
name_V5 = 'V5_mS';% This is a label for the probe to be analysis
bin = 1; % bin of original images, in Hamamatsu sCMOS camera, the bin value determines the average intensity of 
subfolder = 'D:\data\Fused protein expression in stable cell  lines\mS-NES_HEK293T_ASstress\185237_60x confocal\185317_60x confocal\';% a base folder containing the the three tiff stacks (Z-stack) of DAPi, V5 (mS) and G3BP2-EGFP
basepath = '';
pathname = [basepath subfolder];
mkdir([pathname '\analysis']); % saving the preliminary processing results based on the tiff stacks
[file_list] = RScanDir(pathname, '*.tif'); % Using RScanDir function to extract the three tiff stacks (Z-stack) of DAPi, V5 (mS) and G3BP2-EGFP in defined base folder
% load the three tiff file
for n = 1:length(file_list)
    filesize = dir(file_list{n});
    if filesize.bytes <=3.14e+06 %file size who is less than 3 MB needn't be analyzed
        continue
    end
    string_target = 'DAPI';%Keyword 1
    temp = strfind(filesize.name,string_target);
    if isempty('temp') == 0 & temp > 0
        file_DAPI = file_list{n};
    end
    clear temp
    
    string_target = 'V5';%Keyword 2
    temp = strfind(filesize.name,string_target);
    if isempty('temp') == 0 & temp >0
        file_V5= file_list{n};
    end
    clear temp
    
    string_target = 'G3BP2';%Keyword 3
    temp = strfind(filesize.name,string_target);
    if isempty('temp') == 0 & temp >0
        file_G3BP2 = file_list{n};
    end
    clear temp
end

FileTif=[file_DAPI];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov_DAPI=zeros(nImage,mImage,NumberImages,'double');
start_frame = 1;
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov_DAPI(:,:,i)=TifLink.read();
end
TifLink.close();

FileTif=[file_V5];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov_V5=zeros(nImage,mImage,NumberImages,'double');
start_frame = 1;
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov_V5(:,:,i)=TifLink.read();
end
TifLink.close();

FileTif=[file_G3BP2];
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
mov_G3BP2=zeros(nImage,mImage,NumberImages,'double');
start_frame = 1;
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
    TifLink.setDirectory(i);
    mov_G3BP2(:,:,i)=TifLink.read();
end
TifLink.close();
% get the maximum projection of three channels
max_DAPI = double(max(mov_DAPI,[],3));
max_V5 = double(max(mov_V5,[],3));
max_G3BP2 = double(max(mov_G3BP2,[],3));

% RGB converstion of Z-projection images of three channels
saturated = 0.0035;
[cmin_DAPI cmax_DAPI rgb_DAPI]=autorange(max_DAPI,saturated,[0 1 1]);
[cmin_G3BP2 cmax_G3BP2 rgb_G3BP2]=autorange(max_G3BP2,saturated,[1 0 1]);
[cmin_V5 cmax_V5 rgb_V5]=autorange(max_V5,saturated,[0 1 0]);
% display the RGB psuedocolor image of each channel and the overlapped image
figure()
imshow(rgb_DAPI+rgb_G3BP2+rgb_V5,'Border','tight')
saveas(gca,[pathname '\overlay_maxproj_auto.png'])

figure()
imshow(rgb_G3BP2,'Border','tight')
saveas(gca,[pathname '\G3BP2_maxproj_auto.png'])
close(gcf)

figure()
imshow(rgb_V5,'Border','tight')
saveas(gca,[pathname '\V5_maxproj_auto.png'])
close(gcf)

figure()
imshow(rgb_G3BP2,'Border','tight')
% saveas(gca,[pathname '\overlay_maxproj_auto.fig'])
saveas(gca,[pathname '\EGFP_maxproj_auto.png'])
close(gcf)
%% Outline the cells from the merged images of different channels 
% and calculated the average intensity of V5 (mS) and G3BP2-EGFP of the
% selected regions in each slice of the stacks
refimg = 1.5*rgb_DAPI+5*rgb_G3BP2+1.5*rgb_V5;% here we marge the images of 3 channels (DAPI, V5(mS) and G3BP2-EGFP, the factor could be changed for best performance
roi_cells_base = clicky_fullscreenSC(rgb_G3BP2+rgb_G3BP2+rgb_V5, 1, refimg, 0 , 255, 'Outline the cell') % outline one or more cells of interest
saveas(gca,[pathname 'Selected_region.fig'])
saveas(gca,[pathname 'Selected_region.png'])
close(gcf)
roi_bkg = clicky_fullscreenSC(rgb_G3BP2+rgb_G3BP2+rgb_V5, 1, refimg, 0 , 255, 'Outline the background region') % outline a region as background
saveas(gca,[pathname 'Selected_region_bkg.fig'])
saveas(gca,[pathname 'Selected_region_bkg.png'])
close(gcf)
%% calculate the average intensity of V5 (mS) and G3BP2-EGFP of the selected regions in each slice of the stacks
for selected_ROI = 1:length(roi_cells_base) % calculate the regions selected one-by-one
pathname = [subfolder(1:end-1) '_' num2str(selected_ROI)  '\']; 
mkdir([pathname '\analysis']);% save the analysis results of each region to respective folders
roi_cells = roi_cells_base(1,selected_ROI);
[ysize xsize] = size(max_G3BP2);
[x y] = meshgrid(1:xsize, 1:ysize);
inpoly = zeros(ysize,xsize,length(roi_cells));
% the default roi_cells is a 1x1 cell
for n = 1:length(roi_cells)
    inpoly_temp = double(inpolygon(x,y,roi_cells{1,n}(:,1),roi_cells{1,n}(:,2)));% extract selected polygen as a mask
    % calculate the average intensity of G3BP1-EGFP and V5 (mS) from the
    % selected regions slice-by-slice
    intens_V5(:,n) = squeeze(sum(sum(mov_V5.*repmat(inpoly_temp, [1, 1, size(mov_V5, 3)]))))/sum(inpoly_temp(:));
    intens_G3BP2(:,n) = squeeze(sum(sum(mov_G3BP2.*repmat(inpoly_temp, [1, 1, size(mov_G3BP2, 3)]))))/sum(inpoly_temp(:));
    inpoly(:,:,n) = inpoly_temp;
    clear inpoly_temp
end

    % calculate the average background values of G3BP1-EGFP and V5 (mS) from the
    % selected regions slice-by-slice
[ysize xsize] = size(max_G3BP2);
[x y] = meshgrid(1:xsize, 1:ysize);
 
for n = 1:length(roi_cells)
    inpoly_temp = double(inpolygon(x,y,roi_bkg{1,n}(:,1),roi_bkg{1,n}(:,2)));
    intens_V5_bkg(:,n) = squeeze(sum(sum(mov_V5.*repmat(inpoly_temp, [1, 1, size(mov_V5, 3)]))))/sum(inpoly_temp(:));
    intens_G3BP2_bkg(:,n) = squeeze(sum(sum(mov_G3BP2.*repmat(inpoly_temp, [1, 1, size(mov_G3BP2, 3)]))))/sum(inpoly_temp(:));
    inpoly_bkg(:,:,n) = inpoly_temp;
    clear inpoly_temp
end
%% determine the V5-positive pixels
threshold_V5mask_times = 1.5;% the threshold of V5-positve pixels, which is the times of average V5 intensity in each slice 
min_V5_size = 10;% the recognized regions with a area less than is value will be excluded
V5_BW = zeros(ysize,xsize,size(mov_V5,3),length(roi_cells));
pathname_check = [pathname '\check_the_V5_selection\'];
mkdir(pathname_check)
threshold_V5mask = threshold_V5mask_times*(intens_V5-cmin_V5)+cmin_V5;
colororder = interp1([1 2],[400 670],[1:(2-1)/(size(intens_V5,2)-1):2]);
% plot the average intensity of threshold of each slice
figure()
for n = 1:size(intens_V5,2)
    plot(intens_V5(:,n),'-o','Color', getrgb(colororder(n)),'MarkerSize',7);hold on
    plot(threshold_V5mask(:,n),'--o','Color', getrgb(colororder(n)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1 intensity', 'Cell 1 threshold', 'Cell 2 intensity','Cell 2 threshold','Cell 3 intensity','Cell 3 threshold', 'Cell 4 intensity','Cell 4 threshold','Cell 5 intensity','Cell 5 threshold','Cell 6 intensity','Cell 6 threshold','Cell 7 intensity', 'Cell 7 threshold','Cell 8 intensity', 'Cell 8 threshold','Cell 9 intensity', 'Cell 9 threshold','Location','Best')
legend('boxoff')
title(['Average intensity of each slice: RNA channel']);
xlabel('Slice number')
ylabel('Average intensity of each slice')
saveas(gcf,[pathname '\analysis\Cell_intensity_V5.fig'])
saveas(gcf,[pathname '\analysis\Cell_intensity_V5.png'])
close(gcf)

for n = 1:size(mov_V5,3)
    tempslice = mov_V5(:,:,n);
    tempslice = medfilt2(tempslice,[5 5]);%median filter for eliminating salt and pepper noise
    for cell_order = 1:length(roi_cells)
        mask_temp = tempslice.*inpoly(:,:,cell_order);
        BW_temp = mask_temp > threshold_V5mask(n,cell_order);
        
        [L,N] = bwlabel(BW_temp,8);  % 
        % Eliminating the small positive connected components
        for k = 1:N
            puncta_volume(k) = length(find(L == k));
            if puncta_volume(k) <= min_V5_size;
                L(L == k) = 0;
            end
        end
        V5_BW(:,:,n,cell_order) = L;
        %196-213£¬if you run these lines, it will generate the pictures of
        %threshold performance of each slice
%         [cmin_rgbmask,cmax_rgbmask,rgb_mask]=autorange(mask_temp,saturated,[1 1 0]);
%         figure()
%         set(gcf,'outerposition',get(0,'screensize'));
%         subplot(1,3,1)
%         imshow(mask_temp,[cmin_V5 cmax_V5*0.75])
%         title(subplot(1,3,1),['Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(1,3,2)
%         imshow(rgb_mask.*(~BW_temp)+ind2rgb(round(BW_temp*255), pseudocolor(637)))
%         title(subplot(1,3,2),['Check the RNA selection: Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(1,3,3)
%         imshow(rgb_mask.*(~L)+ind2rgb(round(L*255), pseudocolor(637)))
%         title(subplot(1,3,3),['Check the RNA selection after removing small punctas: Cell ' num2str(cell_order) ', slice ' num2str(n) ', threshold = ' num2str(threshold_V5mask_times) 'x avg.'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_V5mask_times) 'x.fig'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_V5mask_times) 'x.png'])
%         close(gcf)
    end
end
V5_BW(V5_BW>0) = 1;% binarize the selcted mask

% generate a Tiff stack for checking the V5-positive mask
for k = 1:size(V5_BW,4)
delete ([pathname_check '\V5_BW_Cell ' num2str(k) '.tif'])
end

for k = 1:size(V5_BW,4)
for n = 1:size(V5_BW,3)
    imwrite(V5_BW(:,:,n,k),[pathname_check '\V5_BW_Cell ' num2str(k) '_1.5.tif'],'WriteMode','append');
end
end
%% determine the G3BP2-positive pixels
threshold_G3BP2mask_times = 5;% ;% the threshold of V5-positve pixels, which is the times of average G3BP2-EGFP intensity in each slice 
min_G3BP2_size = 10;% the recognized regions with a area less than is value will be excluded
G3BP2_BW = zeros(ysize,xsize,size(mov_V5,3),length(roi_cells));
pathname_check = [pathname '\check_the_G3BP2_selection\'];
mkdir(pathname_check)
threshold_G3BP2mask = threshold_G3BP2mask_times*(intens_G3BP2-cmin_G3BP2)+cmin_G3BP2;
colororder = interp1([1 2],[400 670],[1:(2-1)/(size(intens_G3BP2,2)-1):2]);
% plot the average intensity of threshold of each slice
figure()
for n = 1:size(intens_G3BP2,2)
    plot(intens_G3BP2(:,n),'-*','Color', getrgb(colororder(n)),'MarkerSize',7);hold on
    plot(threshold_G3BP2mask(:,n),'--*','Color', getrgb(colororder(n)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1 intensity', 'Cell 1 threshold', 'Cell 2 intensity','Cell 2 threshold','Cell 3 intensity','Cell 3 threshold', 'Cell 4 intensity','Cell 4 threshold','Cell 5 intensity','Cell 5 threshold','Cell 6 intensity','Cell 6 threshold','Cell 7 intensity', 'Cell 7 threshold','Cell 8 intensity', 'Cell 8 threshold','Cell 9 intensity', 'Cell 9 threshold','Location','Best')
legend('boxoff')
title(['Average intensity of each slice: G3BP2 channel']);
xlabel('Slice number')
ylabel('Average intensity of each slice')
saveas(gcf,[pathname '\analysis\Cell_intensity_G3BP2.fig'])
saveas(gcf,[pathname '\analysis\Cell_intensity_G3BP2.png'])
close(gcf)
for n = 1:size(mov_G3BP2,3)
    tempslice = mov_G3BP2(:,:,n);
    tempslice = medfilt2(tempslice,[5 5]);%median filter for eliminating salt and pepper noise
    for cell_order = 1:length(roi_cells)
        mask_temp = tempslice.*inpoly(:,:,cell_order);
        BW_temp = mask_temp > threshold_G3BP2mask(n,cell_order);
        % Eliminating the small positive connected components
 [L,N] = bwlabel(BW_temp,8); 
        for k = 1:N
            puncta_volume(k) = length(find(L == k));
            if puncta_volume(k) <= min_G3BP2_size;
                L(L == k) = 0;
            end
        end
   
        G3BP2_BW(:,:,n,cell_order) = L;
        %269-285£¬if you run these lines, it will generate the pictures of
        %threshold performance of each slice
%         [cmin_rgbmask,cmax_rgbmask,rgb_mask]=autorange(mask_temp,saturated,[1 1 0]);
%         figure()
%         set(gcf,'outerposition',get(0,'screensize'));
%         subplot(1,3,1)
%         imshow(mask_temp,[cmin_G3BP2 cmax_G3BP2])
%         title(subplot(1,3,1),['Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(1,3,2)
%         imshow(rgb_mask.*(~BW_temp)+ind2rgb(round(BW_temp*255), pseudocolor(488)))
%         title(subplot(1,3,2),['Check the G3BP2 selection: Cell ' num2str(cell_order) ', slice ' num2str(n)])
%         subplot(1,3,3)
%         imshow(rgb_mask.*(~L)+ind2rgb(round(L*255), pseudocolor(488)))
%         title(subplot(1,3,3),['Check the G3BP2 selection after removing small punctas: Cell ' num2str(cell_order) ', slice ' num2str(n) ', threshold = ' num2str(threshold_V5mask_times) 'x avg.'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_G3BP2mask_times) 'x.fig'])
%         saveas(gcf,[pathname_check 'Cell ' num2str(cell_order) '_Slice ' num2str(n) '_thresh = ' num2str(threshold_G3BP2mask_times) 'x.png'])
%         close(gcf)
    end
end
G3BP2_BW(G3BP2_BW>0) = 1; % binarize the selcted mask

% generate a Tiff stack for checking the G3BP2-positive mask
for k = 1:size(G3BP2_BW,4)
delete ([pathname_check '\G3BP2_BW_Cell ' num2str(k) '.tif'])
delete ([pathname_check '\Colocalization_BW_Cell ' num2str(k) '.tif'])
end

for k = 1:size(G3BP2_BW,4)
for n = 1:size(G3BP2_BW,3)
    imwrite(G3BP2_BW(:,:,n,k),[pathname_check '\G3BP2_BW_Cell ' num2str(k) '.tif'],'WriteMode','append');
end
end
% generate a Tiff stack for checking the double-positive mask of G3BP2 and V5 (mS) 
for k = 1:size(G3BP2_BW,4)
for n = 1:size(G3BP2_BW,3)
    imwrite(G3BP2_BW(:,:,n,k).*V5_BW(:,:,n,k),[pathname_check '\Colocalization_BW_Cell ' num2str(k) '.tif'],'WriteMode','append');
end
end

%% Calculate the volume (pixel number) of V5 (mS) and G3BP2-EGFP, and the positive ratio of V5 (mS)-positive pixels to G3BP2
for cell_order = 1:size(roi_cells,2)
    inpoly_temp = inpoly(:,:,cell_order);
    cell_size(cell_order) = sum(inpoly_temp(:));
    for n = 1:size(V5_BW,3)       
        V5_positive(cell_order,n) = squeeze(sum(sum(V5_BW(:,:,n,cell_order).*inpoly_temp)));
        V5_in_G3BP2(cell_order,n) = squeeze(sum(sum(V5_BW(:,:,n,cell_order).*inpoly_temp.*G3BP2_BW(:,:,n,cell_order))));
        G3BP2_size(cell_order,n) = sum(sum(G3BP2_BW(:,:,n,cell_order)));
    end
end
%% Calculate and dislplay the volume, average intensity,and total (sum of intesntiy of each postive pixel) of V5 (mS) channel in each slice
[ysize xsize] = size(max_G3BP2);
[x y] = meshgrid(1:xsize, 1:ysize);
for n = 1:length(roi_cells)
   for slice = 1:size(mov_V5,3)
    intens_V5_positive(slice,n) = sum(sum(mov_V5(:,:,slice).*V5_BW(:,:,slice,n)))./sum(sum(V5_BW(:,:,slice,n)))-intens_V5_bkg(slice,n);
    intens_V5_G3BP2(slice,n) = sum(sum(mov_V5(:,:,slice).*V5_BW(:,:,slice,n).*G3BP2_BW(:,:,slice,n)))./sum(sum(V5_BW(:,:,slice,n).*G3BP2_BW(:,:,slice,n)))-intens_V5_bkg(slice,n);
end

end
figure(3)
set(figure(3),'WindowState', 'maximized');
subplot(1,3,1)
for n = 1:size(intens_V5,2)
    plot(V5_positive(n,:),'-o','Color', getrgb(colororder(n)),'MarkerSize',7);hold on
end
box off
axis tight
legend('Cell 1', 'Cell 2','Cell 3', 'Cell 4','Cell 5', 'Cell 6','Cell 7', 'Cell 8','Cell 9', 'Cell 10','Cell 11', 'Cell 12','Cell 13', 'Cell 14','Cell 15', 'Cell 16','Cell 17', 'Cell 18','Cell 19','Location','Best')
legend('boxoff')
title(['Numbers of positive pixel']);
xlabel('Slice number')
ylabel('V5-positive pixel number')
subplot(1,3,2)
plot(intens_V5_positive,'-o','Color', getrgb(colororder(n)),'MarkerSize',7)
box off
axis tight
xlabel('Slice number')
ylabel('intensity (avg.)')
title('V5-positive pixel, avg. inten.')
subplot(1,3,3)
plot(intens_V5_positive.*V5_positive','-o','Color', getrgb(colororder(n)),'MarkerSize',7)
box off
axis tight
xlabel('Slice number')
ylabel('total intensity')
title('V5-positive pixel, total inten.')
saveas(gcf,[pathname '\analysis\V5-positive_distribution.fig'])
saveas(gcf,[pathname '\analysis\V5-positive_distribution.png'])
close(gcf)
%% Select the slices of interest in each cell, we included all slice for calculation
clear cell_size_odd EGFP_positive_odd V5_positive_odd V5_in_G3BP2_odd ratio_in_odd cell_size_even EGFP_positive_even V5_positive_even V5_in_G3BP2_even ratio_in_even
slice_to_be_included = {1:size(intens_G3BP2,1)}; % here we set the default range as all imaged slices

%% Calculate the positive ratio of V5 (mS) to the G3BP2-EGFP in intensity and volume  
% each slice was firstly calculated, and the final results were averaged
% from the previously-selected slices
for n = 1:size(intens_V5,2)
intens_V5_total(n) = sum(intens_V5_positive(slice_to_be_included{n},n).*V5_positive(n,slice_to_be_included{n})','omitnan');
intens_V5_G3BP2_total(n) = sum(intens_V5_G3BP2(slice_to_be_included{n},n).*V5_in_G3BP2(n,slice_to_be_included{n})','omitnan');
volume_V5_G3BP2_selected(n)  = sum(V5_in_G3BP2(n,slice_to_be_included{n}));
volume_V5_total_selected(n) = sum(V5_positive(n,slice_to_be_included{n}));
volume_V5_ratio(n)  = volume_V5_G3BP2_selected(n) ./volume_V5_total_selected(n) ;%Total ratio from the selected slices
intens_V5_ratio(n)  = intens_V5_G3BP2_total(n) ./intens_V5_total(n);%Total ratio from the selected slices
intens_V5_slice_by_slice(:,n) = (intens_V5_G3BP2(:,n).*V5_in_G3BP2(n,:)')./(intens_V5_positive(:,n).*V5_positive(n,:)'); %Ratio in each slice
volume_V5_slice_by_slice(:,n) = (V5_in_G3BP2(n,:)')./(V5_positive(n,:)'); %Ratio in each slice
end
%% Write the results to disk
% write the key results to a xlsx file
filename_to_write = [pathname '\analysis\' name_V5 '_threshold = ' num2str(threshold_V5mask_times) '_G3BP2 threshold = ' num2str(threshold_G3BP2mask_times) '.xlsx'];
delete(filename_to_write)
for n = 1:size(intens_V5_positive,2)
xlswrite(filename_to_write,{'Cell size (max.projection from G3BP2)','Total slice','V5_positive (slice by slice, avg. intensity)','V5_in G3BP2 (slice by slice, avg. intensity)','V5_positive (slice by slice, avg. volume (pixel number))','V5_in G3BP2 (slice by slice,  avg. volume (pixel number))','Selected slices','Sum of V5 on G3BP2 (from selected slices, Volume)','Sum of total (from selected slices, Volume)','Sum of V5 on G3BP2 (from selected slices, Intensity)','Sum of total (from selected slices, Intensity)','Ratio of V5 with G3BP2 (volume)','Ratio of V5 with G3BP2 (Intensity)','threshold of G3BP2','threshold of V5','min_G3BP2_size','min_V5_size'},['Cell ' num2str(n)],['A1']);
xlswrite(filename_to_write,cell_size(n),['Cell ' num2str(n)],['A2']);
xlswrite(filename_to_write,[1:1:size(intens_V5_positive,1)]',['Cell ' num2str(n)],['B2']);
xlswrite(filename_to_write,intens_V5_positive(:,n),['Cell ' num2str(n)],['C2']);
xlswrite(filename_to_write,intens_V5_G3BP2(:,n),['Cell ' num2str(n)],['D2']);
xlswrite(filename_to_write,V5_positive(n,:)',['Cell ' num2str(n)],['E2']);
xlswrite(filename_to_write,V5_in_G3BP2(n,:)',['Cell ' num2str(n)],['F2']);
xlswrite(filename_to_write,slice_to_be_included{1,n}',['Cell ' num2str(n)],['G2']);
xlswrite(filename_to_write,volume_V5_total_selected(n),['Cell ' num2str(n)],['H2']);
xlswrite(filename_to_write,volume_V5_G3BP2_selected(n),['Cell ' num2str(n)],['I2']);
xlswrite(filename_to_write,intens_V5_total(n),['Cell ' num2str(n)],['J2']);
xlswrite(filename_to_write,intens_V5_G3BP2_total(n),['Cell ' num2str(n)],['K2']);
xlswrite(filename_to_write,volume_V5_ratio(n),['Cell ' num2str(n)],['L2']);
xlswrite(filename_to_write,intens_V5_ratio(n),['Cell ' num2str(n)],['M2']);
xlswrite(filename_to_write,threshold_G3BP2mask_times,['Cell ' num2str(n)],['N2']);
xlswrite(filename_to_write,threshold_V5mask_times,['Cell ' num2str(n)],['O2']);
xlswrite(filename_to_write,min_G3BP2_size,['Cell ' num2str(n)],['P2']);
xlswrite(filename_to_write,min_V5_size,['Cell ' num2str(n)],['Q2']);
end

% save the results to a mat file
vars = whos;
saveVars = {vars.name}; % get the variants name list
excludePattern = '^(mov_G3BP2|mov_V5|mov_DAPI)$';% Exclude the variants of imported raw stack data
saveVars = saveVars(cellfun(@isempty, regexp(saveVars, excludePattern))); 
save([pathname 'all variants.mat'],saveVars{:})
end