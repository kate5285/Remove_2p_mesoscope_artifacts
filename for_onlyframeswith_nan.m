%ICholo와 SGholo는 이것으로 run
clear all
mesoh5path = 'D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\h5file_artifactonly\' ;
mesoh5path(mesoh5path=='\') = '/';
mesoSIpath = 'D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\';
mesoSIpath(mesoSIpath=='\') = '/';

foldername =  {'SGholo'};

fnh = [mesoh5path, 'SGholo.h5'];
if exist(fnh, 'file')
    delete(fnh)
end
for f = 1:numel(foldername)
    if f == 1
        tiffns = dir([mesoSIpath, foldername{f}, '/*.tif']);
    else
        tiffns = cat(1, tiffns, dir([mesoSIpath, foldername{f}, '/*.tif']));
    end
end
tiffile = [tiffns(1).folder '\' tiffns(1).name];
tiffheader = imfinfo(tiffile);
hSIh = tiffheader(1).Software;
hSIh = regexp(splitlines(hSIh), ' = ', 'split');
for n=1:length(hSIh)
    if strfind(hSIh{n}{1}, 'SI.hChannels.channelSave')
        nch = n;
        channelssaved = str2num(hSIh{n}{2});
    end
end
numchannels = numel(channelssaved);

artist_info     = tiffheader(1).Artist;
artist_info = artist_info(1:find(artist_info == '}', 1, 'last'));
artist = jsondecode(artist_info);
si_rois = artist.RoiGroups.imagingRoiGroup.rois;
% get ROIs dimensions for each z-plane
nrois = numel(si_rois);

Ly = [];
Lx = [];
cXY = [];
szXY = [];
for k = 1:nrois
    Ly(k,1) = si_rois(k).scanfields(1).pixelResolutionXY(2);
    Lx(k,1) = si_rois(k).scanfields(1).pixelResolutionXY(1);
    cXY(k, [2 1]) = si_rois(k).scanfields(1).centerXY;
    szXY(k, [2 1]) = si_rois(k).scanfields(1).sizeXY;
end
cXY = cXY - szXY/2;
cXY = cXY - min(cXY, [], 1);
mu = median([Ly, Lx]./szXY, 1);
imin = cXY .* mu;

% deduce flyback frames from most filled z-plane
stack = loadFramesBuff(tiffile,1,1,1);

n_rows_sum = sum(Ly);
n_flyback = (size(stack, 1) - n_rows_sum) / max(1, (nrois - 1));

irow = [0 cumsum(Ly'+n_flyback)];
irow(end) = [];
irow(2,:) = irow(1,:) + Ly';

data = struct();
data.nrois = size(irow,2);
for i = 1:size(irow,2)
    data.dx(i) = int32(imin(i,2));
    data.dy(i) = int32(imin(i,1));
    data.lines{i} = irow(1,i):(irow(2,i)-1);
end
numframesperfile = zeros(size(tiffns,1),1);

framecnt = 0;
qcnt = 0;
for f = 1:size(tiffns,1)  
tiffile = [tiffns(f).folder '/' tiffns(f).name];

    reader=ScanImageTiffReader(tiffile);
    desc=reader.descriptions();
    numframesperfile(f) = size(desc,1);

    % temp = reader.data();
    % q = cat(3, q, temp);
    % tic
    q = reader.data();
    edgezonethresh = 1250;
    artifactpixelthresh = 400;
    greenstack=q(:,:,1:2:end);
    redstack = q(:,:,2:2:end);
clearvars q
xedgezone = max(redstack,[],1)>=edgezonethresh; % 1X5166XNframes
% yedgezone = max(redstack,[],2)>=edgezonethresh; % 600X1XNframes
% edgezone = repmat(xedgezone, size(redstack,1),1,1) & repmat(yedgezone, 1,size(redstack,2),1);
% artipixindices = (edgezone==1) & (redstack >= artifactpixelthresh);
artipixindices = (xedgezone == 1);
        for i = 1:size(artipixindices, 2)
            for j = 1:size(artipixindices, 3)
                if artipixindices(1, i, j) == 1
                    artipixindices(1:size(redstack, 1), i, j) = 1;
                end
            end
        end
clearvars edgezone 

[row, col, frame] = ind2sub(size(artipixindices), find(artipixindices == 1));%보정이 필요한 프레임들
clearvars row col
nanframe=unique(frame);
smallredstack = zeros(size(redstack, 1), size(redstack, 2), numel(nanframe));
for i = 1:numel(nanframe)
    frame_idx = nanframe(i);
    smallredstack(:,:,i) = redstack(:,:,frame_idx);
end

smallgreenstack = zeros(size(greenstack, 1), size(greenstack, 2), numel(nanframe));
for i = 1:numel(nanframe)
    frame_idx = nanframe(i);
    smallgreenstack(:,:,i) = greenstack(:,:,frame_idx);
end
%보정이 필요한 프레임 크기로 줄인 본래의 stack들
% clearvars redstack greenstack
newredstack=NaN(size(smallredstack));
xedgezone = max(smallredstack,[],1)>=edgezonethresh; % 이 부분은 일단 작은 스택으로 다시 인덱스 맞추려고 반복
% yedgezone = max(smallredstack,[],2)>=edgezonethresh; 
% edgezone = repmat(xedgezone, size(smallredstack,1),1,1) & repmat(yedgezone, 1,size(smallredstack,2),1);
% artipixindices = (edgezone==1) & (smallredstack >= artifactpixelthresh);%수정이 필요해보임
artipixindices = (xedgezone == 1);
        for i = 1:size(artipixindices, 2)
            for j = 1:size(artipixindices, 3)
                if artipixindices(1, i, j) == 1
                    artipixindices(1:size(redstack, 1), i, j) = 1;
                end
            end
        end
clearvars edgezone 

newredstack(artipixindices==0)=smallredstack(artipixindices==0);%NaN처리, 그 외는 redstack으로 채우기

newgreenstack=NaN(size(smallgreenstack));
newgreenstack(artipixindices==0)=smallgreenstack(artipixindices==0);
clearvars artipixindices
% rimstack9 = NaN(size(smallredstack,1), size(smallredstack,2), size(smallredstack,3), 9);
%         cnt =0;
%         for ii = -1:1
%             for jj = -1:1
%                 cnt = cnt+1;
% 
%                 valx = 2:size(smallredstack,2)-1;
%                 valy = 2:size(smallredstack,1)-1;
%                 rimstack9(valy+ii, valx+jj, :,cnt) = newredstack(valy, valx, :);
%             end
%         end
%         rimstack9=nanmean(rimstack9,4);
%         rimstackinterpol = cat(4, rimstack9);
%         %if there is still nan
%         rimstackinterpol(isnan(rimstackinterpol)) = artifactpixelthresh;
% %green
% gimstack9 = NaN(size(smallgreenstack,1), size(smallgreenstack,2), size(smallgreenstack,3), 9);
%         cnt =0;
%         for ii = -1:1
%             for jj = -1:1
%                 cnt = cnt+1;
% 
%                 valx = 2:size(smallgreenstack,2)-1;
%                 valy = 2:size(smallgreenstack,1)-1;
%                 gimstack9(valy+ii, valx+jj, :,cnt) = newgreenstack(valy, valx, :);
%             end
%         end
%         gimstack9=nanmean(gimstack9,4);
%         gimstackinterpol = cat(4, gimstack9);
%         %if there is still nan
%         gimstackinterpol(isnan(gimstackinterpol)) = artifactpixelthresh;
        interpolstack = NaN(size(newredstack));

        for frame = 1:size(newredstack, 3)
            for i = 1:size(newredstack,1)
            for j = 1:size(newredstack,2)
            if isnan(newredstack(i,j,frame))
            frameIndices = max(1, frame-1):min(size(newredstack, 3), frame+1);
            interpolstack(i, j, frame) = nanmean(newredstack(i, j, frameIndices), 3);
            end
            end
            end
        end
        rimstacktinterpol = newredstack;
        rimstacktinterpol(isnan(newredstack)) = interpolstack(isnan(newredstack));
        for frame = 1:size(newredstack, 3)%NaN이 그래도 남을 경우 한 번 더 interpolation 실행
        for i = 1:size(newredstack,1)
        for j = 1:size(newredstack,2)
        if isnan(rimstacktinterpol(i,j,frame))
            frameIndices = max(1, frame-1):min(size(rimstacktinterpol, 3), frame+1);
            interpolstack(i, j, frame) = nanmean(rimstacktinterpol(i, j, frameIndices), 3);
        end
        end
        end
        end
        rimstacktinterpol(isnan(rimstacktinterpol)) = interpolstack(isnan(rimstacktinterpol));

        interpolstack = NaN(size(newgreenstack));
        for frame = 1:size(newgreenstack, 3)
            for i = 1:size(newgreenstack,1)
            for j = 1:size(newgreenstack,2)
            if isnan(newgreenstack(i,j,frame))
            frameIndices = max(1, frame-1):min(size(newgreenstack, 3), frame+1);
            interpolstack(i, j, frame) = nanmean(newgreenstack(i, j, frameIndices), 3);
            end
            end
            end
        end
        gimstacktinterpol = newgreenstack;
        gimstacktinterpol(isnan(newgreenstack)) = interpolstack(isnan(newgreenstack));
        for frame = 1:size(newgreenstack, 3)%NaN이 그래도 남을 경우 한 번 더 interpolation 실행
        for i = 1:size(newgreenstack,1)
        for j = 1:size(newgreenstack,2)
        if isnan(gimstacktinterpol(i,j,frame))
            frameIndices = max(1, frame-1):min(size(gimstacktinterpol, 3), frame+1);
            interpolstack(i, j, frame) = nanmean(gimstacktinterpol(i, j, frameIndices), 3);
        end
        end
        end
        end
        gimstacktinterpol(isnan(gimstacktinterpol)) = interpolstack(isnan(gimstacktinterpol));
        clearvars interpolstack

 %합치기       
        smallredstack(isnan(newredstack))=rimstacktinterpol(isnan(newredstack));
        smallgreenstack(isnan(newgreenstack))=gimstacktinterpol(isnan(newgreenstack));
%전체에 집어넣기
        for i = 1:numel(nanframe)
        redstack(:,:,nanframe(i)) = smallredstack(:,:,i);
        end
        for i = 1:numel(nanframe)
        greenstack(:,:,nanframe(i)) = smallgreenstack(:,:,i);
        end
        % newstack=zeros(size(redstack,1),size(redstack,2),size(redstack,3)*2);
        clearvars smallredstack smallgreenstack 
        framenumber=size(redstack,3)*2;%초록 빨강 채널 하나로 합치기
        newstack(:,:,2:2:framenumber)=redstack(:,:,:);
        newstack(:,:,1:2:framenumber)=greenstack(:,:,:);
clearvars redstack greenstack

% tempdata = zeros(max(Ly), sum(Lx), size(newstack,3));
    for istrip = 1:numel(Ly)
        tempdata(1:Ly(istrip), sum(Lx(1:istrip-1))+1:sum(Lx(1:istrip)), :) = permute(newstack(:,data.lines{istrip}+1,:),[2 1 3]);
    end
tempdata = permute(tempdata, [2 1 3]);
    k=0;%k설정
% HS: save RESHAPED data in the h5 file!
        qcnt = qcnt+1;
    whodata = whos('tempdata');
    if whodata.bytes <= 4*2^30
        if qcnt ==1
            h5create(fnh,'/data',[size(tempdata,1) size(tempdata,2) Inf], 'DataType','int16', 'ChunkSize',[size(tempdata,1) size(tempdata,2) size(tempdata,3)]); %
            h5write(fnh,'/data',tempdata,[1 1 1],[size(tempdata,1) size(tempdata,2) size(tempdata,3)]);
        sum(isnan(tempdata(:)))
        else
            h5write(fnh,'/data',tempdata,[1 1 framecnt+1],[size(tempdata,1) size(tempdata,2) size(tempdata,3)]);
        sum(isnan(tempdata(:)))
        end
    else
        Nchunks = ceil(whodata.bytes / (4*2^30));
        frameborder = [0:round(size(tempdata,3)/Nchunks):size(tempdata,3) size(tempdata,3)];
        if frameborder(numel(frameborder)) == frameborder(numel(frameborder)-1)
            frameborder(numel(frameborder)) = [];
        end%수정, 정확히 나누어 떨아지는 경우 frameborder에 마지막 두 값이 일치해서 중복되는 거 삭제함.
        for ichunk = 1:length(frameborder)-1
            tempframeinds = frameborder(ichunk)+1:frameborder(ichunk+1); 
            if qcnt == 1
                h5create(fnh,'/data',[size(tempdata,1) size(tempdata,2) Inf], 'DataType','int16', 'ChunkSize',[size(tempdata,1) size(tempdata,2) length(tempframeinds)]);
                h5write(fnh,'/data',tempdata(:,:,tempframeinds),[1 1 1],[size(tempdata,1) size(tempdata,2) length(tempframeinds)]);
            sum(isnan(tempdata(:)))
                qcnt = qcnt+1;%첫번째라서 계속 qcnt==1이라 create랑 write상의 문제가 발생해서 qcnt값을 수정했습니다.
            else
                h5write(fnh,'/data',tempdata(:,:,tempframeinds),[1 1 framecnt+k+1],[size(tempdata,1) size(tempdata,2) length(tempframeinds)]);
            %이 부분에서 오류 찾아서 framecnt이던 것을 length(tempframeinds)로 수정함.
                k=k+length(tempframeinds);
            sum(isnan(tempdata(:))) 
            end
        end
    end
    framecnt = framecnt+size(tempdata,3);
      clearvars newstack tempdata k

    if mod(f,50)==0
    fprintf('%d/%d files loaded\n', f, size(tiffns,1));
    % toc
    end
end
if framecnt ~= sum(numframesperfile)
    error('check frame accumulation')
end
numplanes = 1;
numtimepointsperfile = numframesperfile/(numchannels*numplanes);

save(strcat(mesoh5path, 'h5params.mat'), 'foldername', 'tiffns', ...
    'numframesperfile', 'numtimepointsperfile', 'numchannels') %, 'whichvisfile'
disp(foldername)
disp(framecnt)

%%video
file_path ="D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\h5file_artifactonly\ICholo.h5";  % HDF5 파일 경로 설정
dataset_name = '/data';  % 데이터셋 경로 설정
frame_range = 12500:13000;
image_data = h5read(file_path, dataset_name, [1, 1, frame_range(1)], [Inf, Inf,numel(frame_range)]);%Sgholo같이 너무 많은 경우만
import matlab.io.VideoWriter.*
% rolling average video,red one

adjustedStack = image_data;

adjustedStack(adjustedStack <= -20) = -20;
adjustedStack(adjustedStack >= 150) = 150;
adjustedStack = (adjustedStack + 20) / 170;
adjustedStack=permute(adjustedStack(:,:,1:2:end), [2 1 3]);
rollframes = 1; 
rollind = (1:rollframes)'+[0:size(adjustedStack,3)-rollframes];
vr=VideoWriter('sgholo_lines_fov_green.mp4', 'MPEG-4');
vr.FrameRate=18.2;
caxis ([0,150]);
open(vr);
for ii = 1:size(rollind,2)
writeVideo(vr, squeeze(mean(adjustedStack(:,:,rollind(:,ii)), 3) ));
end
close(vr);