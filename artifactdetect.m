import ScanImageTiffReader.ScanImageTiffReader.*
tiffDirectory = 'D:\doyeon_kim\MesoHoloExpts_scanimage_MU31_2_ICholo\230106';
tiffFiles = dir(fullfile(tiffDirectory, '*.tif'));

redthresh = 512;
nTrials = 50;%numel(tiffFiles);

imstack = [];
for trial = 1:nTrials
    tiffn = fullfile(tiffDirectory, tiffFiles(trial).name);
    imfile = ScanImageTiffReader(tiffn);
    imdata = imfile.data();

    nFrames = size(imdata, 3);
    imstack = cat(3, imstack, imdata);
end
imstack = double(imstack);

greenstack = imstack(:,:,1:2:end);
redstack = imstack(:,:,2:2:end);

redstack(redstack>redthresh) = NaN;
greenstack (greenstack>redthresh) = NaN;
% edgezone
xedgezone = max(redstack,[],1); % 1X5166XNframes
yedgezone = max(redstack,[],2); % 600X1XNframes
edgezone = repmat(xedgezone, size(redstack,1),1,1) & repmat(yedgezone, 1,size(redstack,2),1);

exredframe = squeeze(redstack(:,:,8));
exedgezone = max(exredframe,[],1)>512 & max(exredframe,[],2)>512;
figure; imshow(exedgezone)
%replace NaN values_spatially_Replace the NaN values as the average of the pixels around it in space
%for rimstack9 method

% rimstack9 = NaN(size(redstack,1), size(redstack,2), size(redstack,3), 9);
% cnt =0;
% for ii = -1:1
% for jj = -1:1
% cnt = cnt+1;
% rimstack9 = NaN(size(redstack));
% valx = 2:size(redstack,2)-1;
% valy = 2:size(redstack,1)-1;
% rimstack9(valy+ii, valx+jj, :,cnt) = redstack(valy, valx, :);
% end
% end
% rimstack9=nanmean(rimstack9,4);
% rimstackinterpol = cat(4, rimstack9);
% rimstacktinterpol(isnan(rimstacktinterpol)) = 512;


avg_vals = NaN(size(redstack));
for frame = 1:size(redstack, 3)
    for ii = 1:size(redstack, 1)
        for jj = 1:size(redstack, 2)
             
            % NaN인 픽셀에 대해서만 평균 계산
         if isnan(redstack(ii, jj,frame))
            
            if ii == 1 && size(redstack, 2)>jj &&jj>1
                neighbors=[redstack(ii, jj-1, frame), redstack(ii, jj+1, frame), ...
                           redstack(ii+1, jj-1, frame), redstack(ii+1, jj, frame), redstack(ii+1, jj+1, frame)];
            
            elseif jj==1 && size(redstack, 1)>ii&& ii>1
                neighbors=[redstack(ii-1, jj, frame), redstack(ii-1, jj+1, frame), ...
                           redstack(ii, jj+1, frame),redstack(ii+1, jj, frame), redstack(ii+1, jj+1, frame)];
            
            elseif ii==1 && jj==1
                neighbors=[redstack(ii, jj+1, frame),redstack(ii+1, jj, frame),redstack(ii+1, jj+1, frame)];
            
            elseif ii == 1 && size(redstack, 2)==jj
                neighbors=[redstack(ii, jj-1, frame), ...
                           redstack(ii+1, jj-1, frame), redstack(ii+1, jj, frame)];
            
            elseif jj==size(redstack, 2) && size(redstack, 1)>ii&& ii>1
                neighbors=[redstack(ii-1, jj, frame), redstack(ii-1, jj-1, frame),redstack(ii, jj-1, frame) ...
                           redstack(ii+1, jj, frame), redstack(ii+1, jj-1, frame)];
            
            elseif ii==size(redstack, 1) && jj==1
                neighbors=[redstack(ii, jj+1, frame),redstack(ii-1, jj, frame),redstack(ii-1, jj+1, frame)];
            
            elseif ii==size(redstack, 1) && jj==size(redstack, 2)
                neighbors=[redstack(ii, jj-1, frame),redstack(ii-1, jj, frame),redstack(ii-1, jj-1, frame)];
            
            elseif ii == size(redstack, 1) && size(redstack, 2)>jj &&jj>1
                neighbors=[redstack(ii, jj-1, frame), redstack(ii, jj+1, frame), ...
                           redstack(ii-1, jj-1, frame), redstack(ii-1, jj, frame), redstack(ii-1, jj+1, frame)];
            
            elseif size(redstack, 1)>ii&&ii>1 && size(redstack, 2)>jj&&jj>1
                % 인접한 8개 픽셀의 값들을 가져옴
                neighbors = [redstack(ii-1, jj-1, frame), redstack(ii-1, jj, frame), redstack(ii-1, jj+1, frame), ...
                             redstack(ii, jj-1, frame), redstack(ii, jj+1, frame), ...
                             redstack(ii+1, jj-1, frame), redstack(ii+1, jj, frame), redstack(ii+1, jj+1, frame)];
            end
                % NaN이 아닌 값들로 평균 계산
                avg = nanmean(neighbors);
                avg_vals(ii, jj,frame) = avg; 
         end    
        end
        
    end
end
rimstackinterpol = redstack;
rimstackinterpol(isnan(redstack)) = avg_vals(isnan(redstack));
%green spatial
avg_vals = NaN(size(greenstack));
for frame = 1:size(greenstack, 3)
    for ii = 1:size(greenstack, 1)
        for jj = 1:size(greenstack, 2)
             
            % NaN인 픽셀에 대해서만 평균 계산
         if isnan(greenstack(ii, jj,frame))
            
            if ii == 1 && size(greenstack, 2)>jj &&jj>1
                neighbors=[greenstack(ii, jj-1, frame), greenstack(ii, jj+1, frame), ...
                           greenstack(ii+1, jj-1, frame), greenstack(ii+1, jj, frame), greenstack(ii+1, jj+1, frame)];
            
            elseif jj==1 && size(greenstack, 1)>ii&& ii>1
                neighbors=[greenstack(ii-1, jj, frame), greenstack(ii-1, jj+1, frame), ...
                           greenstack(ii, jj+1, frame),greenstack(ii+1, jj, frame), greenstack(ii+1, jj+1, frame)];
            
            elseif ii==1 && jj==1
                neighbors=[greenstack(ii, jj+1, frame),greenstack(ii+1, jj, frame),greenstack(ii+1, jj+1, frame)];
            
            elseif ii == 1 && size(greenstack, 2)==jj
                neighbors=[greenstack(ii, jj-1, frame), ...
                           greenstack(ii+1, jj-1, frame), greenstack(ii+1, jj, frame)];
            
            elseif jj==size(greenstack, 2) && size(greenstack, 1)>ii&& ii>1
                neighbors=[greenstack(ii-1, jj, frame), greenstack(ii-1, jj-1, frame),greenstack(ii, jj-1, frame) ...
                           greenstack(ii+1, jj, frame), greenstack(ii+1, jj-1, frame)];
            
            elseif ii==size(greenstack, 1) && jj==1
                neighbors=[greenstack(ii, jj+1, frame),greenstack(ii-1, jj, frame),greenstack(ii-1, jj+1, frame)];
            
            elseif ii==size(greenstack, 1) && jj==size(greenstack, 2)
                neighbors=[greenstack(ii, jj-1, frame),greenstack(ii-1, jj, frame),greenstack(ii-1, jj-1, frame)];
            
            elseif ii == size(greenstack, 1) && size(greenstack, 2)>jj &&jj>1
                neighbors=[greenstack(ii, jj-1, frame), greenstack(ii, jj+1, frame), ...
                           greenstack(ii-1, jj-1, frame), greenstack(ii-1, jj, frame), greenstack(ii-1, jj+1, frame)];
            
            elseif size(greenstack, 1)>ii&&ii>1 && size(greenstack, 2)>jj&&jj>1
                % 인접한 8개 픽셀의 값들을 가져옴
                neighbors = [greenstack(ii-1, jj-1, frame), greenstack(ii-1, jj, frame), greenstack(ii-1, jj+1, frame), ...
                             greenstack(ii, jj-1, frame), greenstack(ii, jj+1, frame), ...
                             greenstack(ii+1, jj-1, frame), greenstack(ii+1, jj, frame), greenstack(ii+1, jj+1, frame)];
            end
                % NaN이 아닌 값들로 평균 계산
                avg = nanmean(neighbors);
                avg_vals(ii, jj,frame) = avg; 
         end    
        end
        
    end
end
gimstackinterpol = greenstack;
gimstackinterpol(isnan(greenstack)) = avg_vals(isnan(greenstack));
%need to check why green has nan. 주위 nan이 많은 경우인듯
gimstackinterpol(isnan(gimstackinterpol)) = redthresh;

%for temporal replacement red
nans = NaN(size(redstack));
nanMask = isnan(redstack);
for frame = 1:size(redstack, 3)
    for i = 1:size(redstack,1)
    for j = 1:size(redstack,2)
    if nanMask (i,j,frame)
    frameIndices = max(1, frame-1):min(size(redstack, 3), frame+1);
    nans(i, j, frame) = nanmean(redstack(i, j, frameIndices), 3);
    end
    end
    end
end
rimstacktinterpol = redstack;
rimstacktinterpol(isnan(redstack)) = nans(isnan(redstack));
rimstacktinterpol(isnan(rimstacktinterpol)) = redthresh;

%red y coordinates PMT artifact extract
% x_range = [1:140, 460:600];
% nanMask = isnan(redstack);
% nanycoords = cell(size(redstack, 3), length(x_range));
% for frame = 1:size(nanMask, 3)
% for ix = 1:length(x_range)
%         x = x_range(ix);
%         nanYCoords = find(nanMask(x, :, frame));
%     nanycoords{frame, ix} = nanYCoords;
% 
% end
% end
% 
% nanp = NaN(size(redstack));
% 
% for frame = 1:size(redstack, 3)
%     for i = 140:460
%         nanYCoordsFrame = nanycoords{frame, :};
%         for j = nanYCoordsFrame
%     frameIndices = max(1, frame-3):min(size(redstack, 3), frame+3);
%     nanp(i, j, frame) = nanmean(redstack(i, j, frameIndices));
%         end
%     end
% end
% rpmttinterpol = nanp;
% rpmttinterpol(isnan(nanp)) = redstack(isnan(nanp));
% rpmttinterpol(isnan(rpmttinterpol)) = redthresh;

%방안 2

% redthresh=1000;
% redstack(redstack>redthresh) = NaN;
redthresh1=550;
nanMask = isnan(redstack);
nank = NaN(size(redstack));
nanycoords = cell(size(redstack, 3), size(redstack, 1));
for frame = 1:size(nanMask, 3)
for x = 1:size(nanMask,1)
    nanYCoords = find(nanMask(x, :, frame));
    nanycoords{frame, x} = nanYCoords;
end
end
%using 300indexes-> FIXED TO ENTIRE X INDEX
for frame = 1:size(redstack, 3)
    for i = 1:size(redstack, 1)
        nanYCoordsFrame = horzcat(nanycoords{frame, :});
        nanYCoordsFrame=unique(nanYCoordsFrame);
        for k = 1: floor(numel(nanYCoordsFrame) / 300)
            startIdx = 1 + (k - 1) * 300;
            endIdx = min(nanYCoordsFrame(startIdx) + 299,nanYCoordsFrame(end));
            nanYCoordsSubset = nanYCoordsFrame(startIdx):endIdx;
        for  j = nanYCoordsSubset
    frameIndices = max(1, frame+7):min(size(redstack, 3), frame+10);
    nank(i, j, frame) = nanmean(redstack(i, j, frameIndices));
        end
        end
    end
end
rpmttinterpol = redstack;
rpmttinterpol(isnan(redstack)) = nank(isnan(redstack));
%for those not interpolated
rpmttinterpol(isnan(rpmttinterpol)) = redthresh1;

%what if all nanycoorframes
for frame = 1:size(redstack, 3)
    for i = 1:size(redstack, 1)
        nanYCoordsFrame = horzcat(nanycoords{frame, :});
        nanYCoordsFrame=unique(nanYCoordsFrame);
        for  j = nanYCoordsFrame
    frameIndices = max(1, frame+7):min(size(redstack, 3), frame+10);
    nank(i, j, frame) = nanmean(redstack(i, j, frameIndices));
        end
    end
end
rpmttinterpol = nank;
rpmttinterpol(isnan(nank)) = redstack(isnan(nank));
%for those not interpolated
rpmttinterpol(isnan(rpmttinterpol)) = redthresh1;

%what if start to end ycoords,, I don't think this is a good idea
% for frame = 1:size(redstack, 3)
%     for i = 140:460
%         nanYCoordsFrame = horzcat(nanycoords{frame, :});
%         nanYCoordsFrame=unique(nanYCoordsFrame);
%         startidx=round(nanYCoordsFrame(1));
%         endidx=round(nanYCoordsFrame(end));
%         for  j = startidx:endidx
%     frameIndices = max(1, frame+6):min(size(redstack, 3), frame+10);
%     nank(i, j, frame) = nanmean(redstack(i, j, frameIndices));
%         end
%     end
% end
% rpmttinterpol = nank;
% rpmttinterpol(isnan(nank)) = redstack(isnan(nank));
%green
nans = NaN(size(greenstack));
nanMask = isnan(greenstack);
for frame = 1:size(greenstack, 3)
    for i = 1:size(greenstack,1)
    for j = 1:size(greenstack,2)
    if nanMask (i,j,frame)
    frameIndices = max(1, frame-1):min(size(greenstack, 3), frame+1);
    nans(i, j, frame) = nanmean(greenstack(i, j, frameIndices), 3);
    end
    end
    end
end
gimstacktinterpol = greenstack;
gimstacktinterpol(isnan(greenstack)) = nans(isnan(greenstack));
gimstacktinterpol(isnan(gimstacktinterpol)) = redthresh;

% -20 ~ 150 까지로 제한. 150 이상은 150으로 대체, -20 이하는 -20으로 대체, 0 - 1 range로 normalize 
%do it for rimstackinterpol and rimstacktinterpol
Nframes=size(redstack,3);
resizedgreenstack = zeros(1200,2400,1255);
resizedredstack = zeros(1200,2400,1255);
% for istrip = 1:4
% resizedgreenstack (1:600,1200*(istrip-1)+1:1200*(istrip),:)=gimstackinterpol(:,(istrip-1)*1322+1:(istrip)*1200+122*(istrip-1),:);
% resizedredstack (1:600,1200*(istrip-1)+1:1200*(istrip),:)=rimstackinterpol(:,(istrip-1)*1322+1:(istrip)*1200+122*(istrip-1),:);
% end
for istrip = 1:4
resizedgreenstack(1:1200,600*(istrip-1)+1:600*(istrip),:)=permute(gimstackinterpol(:, (istrip-1)*1322+1:(istrip)*1200+122*(istrip-1), :), [2 1 3]);
resizedredstack(1:1200,600*(istrip-1)+1:600*(istrip),:)=permute(rimstackinterpol(:, (istrip-1)*1322+1:(istrip)*1200+122*(istrip-1), :), [2 1 3]);
end
%make video
import matlab.io.VideoWriter.*
% rolling average video,red one

adjustedStack = resizedredstack;

adjustedStack(adjustedStack <= -20) = -20;
adjustedStack(adjustedStack >= 150) = 150;
adjustedStack = (adjustedStack + 20) / 170;

rollframes = 1; 
rollind = (1:rollframes)'+[0:size(adjustedStack,3)-rollframes];
vr=VideoWriter('trials50redroll1raw(1)_fov.mp4', 'MPEG-4');
vr.FrameRate=18.2;
caxis ([0,150]);
open(vr);
for ii = 1:size(rollind,2)
writeVideo(vr, squeeze(mean(adjustedStack(:,:,rollind(:,ii)), 3) ));
end
close(vr);

%green one
adjustedStack = resizedgreenstack;
adjustedStack(adjustedStack <= -20) = -20;
adjustedStack(adjustedStack >= 150) = 150;
adjustedStack = (adjustedStack + 20) / 170;

rollframes = 4; 
rollind = (1:rollframes)'+[0:size(adjustedStack,3)-rollframes];
vg=VideoWriter('trials50greenroll4temp_fov.mp4', 'MPEG-4');
vg.FrameRate=18.2;
caxis([0, 150]);
open(vg);
for ii = 1:size(rollind,2)
writeVideo(vg, squeeze( mean(adjustedStack(:,:,rollind(:,ii)), 3) ));
end
close(vg);

%to erase artifact on edges after erasing them in the middle
% redthresh1=550;
rmptsides=rpmttinterpol;
rmptsides(rmptsides>=redthresh1) = NaN;

x_range = [1:160, 440:600];
nans = NaN(size(rmptsides));
nanMaskt = isnan(rmptsides);
for frame = 1:size(rmptsides, 3)
    for i = x_range
    for j = 1:size(rmptsides,2)
    if nanMaskt (i,j,frame)
    frameIndices = max(1, frame+7):min(size(rmptsides, 3), frame+9);
    nans(i, j, frame) = nanmean(rmptsides(i, j, frameIndices));
    end
    end
    end
end
rimstacktinterpol(x_range,:,:) =nans(x_range,:,:);
rimstacktinterpol(160:440,:,:) = nank(160:440,:,:);
rimstacktinterpol(isnan(rimstacktinterpol))=redstack(isnan(rimstacktinterpol));
rimstacktinterpol(isnan(rimstacktinterpol)) = redthresh1;
