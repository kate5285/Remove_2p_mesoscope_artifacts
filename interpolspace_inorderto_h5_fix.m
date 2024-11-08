% foldername = {'ICholo','SGholo','stimtest_5cph'};%전체 폴더
% for folder= 1: numel(foldername)
% clearvars -except folder;
clear all
mesoh5path = 'D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\h5file_artifactonly\' ;
mesoh5path(mesoh5path=='\') = '/';
mesoSIpath = 'D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\';
mesoSIpath(mesoSIpath=='\') = '/';
ls(mesoSIpath)
%%
sireaderpath='d:\Users\USER\Documents\MATLAB';
addpath(genpath(sireaderpath))
import ScanImageTiffReader.ScanImageTiffReader.*;

foldername = {'SGholo'};%전체 폴더

% Name='file00011';
folder=1;%nedd to change
fnh = [mesoh5path, foldername{folder}, '.h5'];
if exist(fnh, 'file')
    delete(fnh)
end
% HS: add tif files that are not in a folder
% tiffns = dir([mesoSIpath,  '/*.tif']);
% for f = 1:numel(foldername)
f=1;%need to chaange
    if f == 1
        tiffns = dir([mesoSIpath, foldername{folder}, '/*.tif']);
    else
        tiffns = cat(1, tiffns, dir([mesoSIpath, foldername{folder}, '/*.tif']));
    end
% end

%% dimensions of mesoscope scanimage tif files
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

%% load all tif files, reshape, crop and save as h5
% with Nframes2read = 10000, every qcnt takes 1.5 min
% still estimated ~45min for a holo session
% when calling h5write for every tif file, this takes a *long* time.
% approx 10 files per minute, when each file contains 200 frames (10000 frames ~5 min)
% tic
numframesperfile = zeros(size(tiffns,1),1);

framecnt = 0;
qcnt = 0;

% whichvisfile = zeros(size(suite2pfilelist,1),1);
for f = 1:size(tiffns,1)   
    tiffile = [tiffns(f).folder '/' tiffns(f).name];

    reader=ScanImageTiffReader(tiffile);
    desc=reader.descriptions();
    numframesperfile(f) = size(desc,1);

    % temp = reader.data();
    % q = cat(3, q, temp);
    % tic
    q = reader.data();

    temptiffolder = tiffns(f).folder;
    temptiffolder(tiffns(f).folder=='\') = '/';
    tempfolders = strsplit(temptiffolder, mesoSIpath);
    if numel(tempfolders)==1
        tiffolder = {};
    else
        tiffolder = tempfolders{2};
    end
    % HS: remove artifact only for folders that contain 'holo' or 'stimtest'
    if any (contains(tiffolder, 'holo')) || any (contains(tiffolder, 'stimtest')) %배열 문제 때문에 any를 붙임.
        edgezonethresh = 1250;%원래 512
        % artifactpixelthresh = 400;%원래 400
        greenstack=q(:,:,1:2:end);
        redstack = q(:,:,2:2:end);
        clearvars q %대용량 파일의 경우만
        % red channel value가 edgezone 이면서 artifactpixelthresh 넘어가는 픽셀들을 NaN으로 처리.
        % -> artifact pixel들이 nan으로 되어있는 stack을 newgreenstack, newredstack이라 부르자
        xedgezone = max(redstack,[],1)>=edgezonethresh; % 1X5166XNframes
        % yedgezone = max(redstack,[],2)>=edgezonethresh; % 600X1XNframes
        % edgezone = repmat(xedgezone, size(redstack,1),1,1) & repmat(yedgezone, 1,size(redstack,2),1);

        % artipixindices = (edgezone==1) & (redstack >= artifactpixelthresh);
        %artipixindices를 사용할거면        
        artipixindices = (xedgezone == 1);
        for i = 1:size(artipixindices, 2)
            for j = 1:size(artipixindices, 3)
                if artipixindices(1, i, j) == 1
                    artipixindices(1:size(redstack, 1), i, j) = 1;
                end
            end
        end
        %그 후에
        newredstack=NaN(size(redstack));
        newredstack(artipixindices==0)=redstack(artipixindices==0);
        clearvars redstack
        newgreenstack=NaN(size(greenstack));
        newgreenstack(artipixindices==0)=greenstack(artipixindices==0);
        clearvars greenstack artipixindices
        % artipixindices = (xedgezone==1);
        % [rowIndices, columnIndices, frameIndices] = ind2sub(size(artipixindices), find(artipixindices));
        % coordinates= [rowIndices, columnIndices, frameIndices];
        % newredstack = double(redstack);
        % for i = 1:size(coordinates, 1)
        %     % row= coordinates(i,1);
        %     col = coordinates(i,2);
        %     frame = coordinates(i,3);
        %     % newredstack(row, col, frame) = NaN;
        %     newredstack(:, col, frame) = NaN;
        % end
        %green
        % newgreenstack = double(greenstack);
        % for i = 1:size(coordinates, 1)
        %    % row= coordinates(i,1);
        %     col = coordinates(i,2);
        %     frame = coordinates(i,3);
        %     newgreenstack(:, col, frame) = NaN;
        % end

        %interpol for red- spatial 
        % rimstack9 = NaN(size(redstack,1), size(redstack,2), size(redstack,3), 9);
        % cnt =0;
        % for ii = -1:1
        %     for jj = -1:1
        %         cnt = cnt+1;
        % 
        %         valx = 2:size(redstack,2)-1;
        %         valy = 2:size(redstack,1)-1;
        %         rimstack9(valy+ii, valx+jj, :,cnt) = newredstack(valy, valx, :);
        %     end
        % end
        % rimstack9=nanmean(rimstack9,4);
        % rimstackinterpol = cat(4, rimstack9);
        % %if there is still nan
        % rimstackinterpol(isnan(rimstackinterpol)) = artifactpixelthresh;
        % 
        % %interpol for green
        % gimstack9 = NaN(size(greenstack,1), size(greenstack,2), size(greenstack,3), 9);
        % cnt =0;
        % for ii = -1:1
        %     for jj = -1:1
        %         cnt = cnt+1;
        % 
        %         valx = 2:size(greenstack,2)-1;
        %         valy = 2:size(greenstack,1)-1;
        %         gimstack9(valy+ii, valx+jj, :,cnt) = newgreenstack(valy, valx, :);
        %     end
        % end
        % gimstack9=nanmean(gimstack9,4);
        % gimstackinterpol = cat(4, gimstack9);
        % %if there is still nan
        % gimstackinterpol(isnan(gimstackinterpol)) = artifactpixelthresh;

        %for temporal replacement red
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
        % newgreenstack, newredstack에서 NaN으로 되어있는 픽셀들에 한해서 동일 위치의
        % rimstackinterpol, gimstackinterpol 값으로 대체
        % red=q(:,:,2:2:end);
        newredstack(isnan(newredstack))=rimstacktinterpol(isnan(newredstack));
        clearvars rimstacktinterpol
        % green=q(:,:,1:2:end);
        newgreenstack(isnan(newgreenstack))=gimstacktinterpol(isnan(newgreenstack));
        clearvars gimstacktinterpol
        newstack=zeros(size(newredstack,1),size(newredstack,2),size(newredstack,3)*2);
        newstack(:,:,2:2:end)=newredstack;
        newstack(:,:,1:2:end)=newgreenstack;
    % else
    %     newstack = q;
    end

    % 여기서 h5에 덧붙여 저장하기
    % reshape the data
    tempdata = zeros(max(Ly), sum(Lx), size(newstack,3));
    for istrip = 1:numel(Ly)
        tempdata(1:Ly(istrip), sum(Lx(1:istrip-1))+1:sum(Lx(1:istrip)), :) = permute(newstack(:,data.lines{istrip}+1,:),[2 1 3]);
    end
        clearvars newstack
tempdata = permute(tempdata, [2 1 3]);    %여기가 안 된 거 같음
    k=0;%k설정
% HS: save RESHAPED data in the h5 file!
    qcnt = qcnt+1;
    whodata = whos('tempdata');
    if whodata.bytes <= 4*2^30
        if qcnt ==1
            h5create(fnh,'/data',[size(tempdata,1) size(tempdata,2) Inf], 'DataType','int16', 'ChunkSize',[size(tempdata,1) size(tempdata,2) size(tempdata,3)]); %
            h5write(fnh,'/data',tempdata,[1 1 1],[size(tempdata,1) size(tempdata,2) size(tempdata,3)]);
        sum(isnan(tempdata(:)))%혹시 NAN이 있는지
        else
            h5write(fnh,'/data',tempdata,[1 1 framecnt+1],[size(tempdata,1) size(tempdata,2) size(tempdata,3)]);
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
                h5write(fnh,'/data',tempdata(:,:,tempframeinds+k),[1 1 1],[size(tempdata,1) size(tempdata,2) length(tempframeinds)]);
            sum(isnan(tempdata(:)))
                qcnt = qcnt+1;%첫번째라서 계속 qcnt==1이라 create랑 write상의 문제가 발생해서 qcnt값을 수정했습니다.
            k=k+length(tempframeinds);
            else
                h5write(fnh,'/data',tempdata(:,:,tempframeinds),[1 1 framecnt+k+1],[size(tempdata,1) size(tempdata,2) length(tempframeinds)]);
            %이 부분에서 오류 찾아서 framecnt이던 것을 length(tempframeinds)로 수정함.
                k=k+length(tempframeinds);
             end
        end
    end
    framecnt = framecnt+size(tempdata,3);

    %newstack, q, tempdata 비우기
    clearvars tempdata k

    if mod(f,50)==0
        fprintf('%d/%d files loaded\n', f, size(tiffns,1));
        % toc
    end
end
% end
% toc
if framecnt ~= sum(numframesperfile)
    error('check frame accumulation')
end
numplanes = 1;
numtimepointsperfile = numframesperfile/(numchannels*numplanes);

save(strcat(mesoh5path, 'h5params.mat'), 'foldername', 'tiffns', ...
    'numframesperfile', 'numtimepointsperfile', 'numchannels') %, 'whichvisfile'
disp(Name)
disp(framecnt)
% end