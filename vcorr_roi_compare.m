clear all 
close all
mousedate = 'MU31_2/230106/';
drivepath = '\\shinlab\ShinLab\MesoHoloExpts\';
mesoSIpath = [drivepath 'mesoholoexpts_scanimage/' mousedate];
onlinepath = [drivepath 'mesoholoexpts_scanimage/' mousedate 'ClosedLoop_justgreen/'];
newpath = '\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts\MU31_2\230106\suite2p\combined\';
%newpath = 'D:\doyeon_kim\previous_data_suite2p\';
pathpp = [drivepath 'mesoholoexpts_postprocessed/' mousedate];

news2p = load([newpath 'Fall.mat']);

cd([mesoSIpath 'visstiminfo'])
V = dir([mesoSIpath 'visstiminfo/*.mat']);
Vname = {V.name};
orderv = zeros(numel(V),1);
for ii = 1:numel(V)
    tempVsplit = strsplit(Vname{ii}, '.mat');
    orderv(ii) = str2double(tempVsplit{1}(end-2:end));
end
[~,si] = sort(orderv);
Vordered = Vname(si);
disp(Vordered')

dircontents = dir(mesoSIpath);
mesoSIfolders = {dircontents([dircontents.isdir]).name};
disp(mesoSIfolders')
open mesoSIfolders

compute_numframes = true;
allfoldernames = {'101','retinotopy0','staticICtxi0','staticgratings','staticgratings12', ...
    'staticICtxi1', 'RFcircleC', 'sizecircleC', 'stimtest_5cph', 'SGholo', 'ICholo'};
nexpts = {'x', '000', '001', '002', '003', '004', '005', '006', 'x', 'x', 'x'};

if ~( all(ismember(allfoldernames, {dircontents.name})) )
    disp(allfoldernames(~ismember(allfoldernames, {dircontents.name})) )
    error('check foldername: user did not input correctly')
end

if ~exist(pathpp, 'dir')
    mkdir(pathpp)
end
%%
[holoreqfile,holoreqpath] = uigetfile([onlinepath '*.mat'], 'holoRequest_MU31_2_230106_staticICtxi_001');
load([holoreqpath holoreqfile])
fullnpix = [sum(nxpix),mean(nypix)];
[fullxpix,fullypix] = meshgrid(linspace(1,fullnpix(1),fullnpix(1)),...
    linspace(1,fullnpix(2),fullnpix(2)));
fullxsize = sum(xsize); %um
fullysize = mean(ysize); %um
fullxcenter = (xcenter(1)-xsize(1)/2) + ...
    abs((xcenter(1)-xsize(1)/2) - (xcenter(end)+xsize(end)/2))/2; %um
fullycenter = mean(ycenter); %um
xumperpix = fullxsize/size(news2p.ops.meanImg,2);
yumperpix = fullysize/size(news2p.ops.meanImg,1);
%%
cropedgethr = 0;
[reftif, refpath] = uigetfile('\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts_scanimage\MU31_2\230106\sizecircleC\file_00210.tif');
fname = [refpath reftif];
header = imfinfo(fname);
artist_info     = header(1).Artist;
artist_info = artist_info(1:find(artist_info == '}', 1, 'last'));
artist = jsondecode(artist_info);
hSIh = header(1).Software;
hSIh = regexp(splitlines(hSIh), ' = ', 'split');
for n=1:length(hSIh)
	if strfind(hSIh{n}{1}, 'SI.hRoiManager.scanVolumeRate')
		fs = str2double(hSIh{n}{2});
	end
	if strfind(hSIh{n}{1}, 'SI.hFastZ.userZs')
		zs = str2num(hSIh{n}{2}(2:end-1));
		nplanes = numel(zs);
	end
	
end
si_rois = artist.RoiGroups.imagingRoiGroup.rois;
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

if ~(sum(Lx)==news2p.ops.Lx && isequal(unique(Ly), news2p.ops.Ly))
    error('unexpected Lx/Ly error')
end
%%
load([onlinepath 'online_params.mat'])
online = load([onlinepath '/suite2p/plane0/Fall.mat']);
%vcorr 떨어진 거 확인용
newyfar=(1200-size(news2p.ops.Vcorr,1))/2;%1200은 meanimg의 y길이
newxfar=(2400-size(news2p.ops.Vcorr,2))/2;
news2pvcorr=zeros(news2p.ops.Ly, news2p.ops.Lx);
news2pvcorr(newyfar+1:news2p.ops.Ly-newyfar, newxfar+1:news2p.ops.Lx-newxfar)=news2p.ops.Vcorr;
newstat = cell2mat(news2p.stat);
newmed = double( cat(1, newstat.med) );
%% 'iscell' criterion for pre-manual curation
tic
% mean Vcorr per ROIs
Nnewrois = numel(news2p.stat);
imnewroi=zeros(news2p.ops.Ly, news2p.ops.Lx);
newroictr = zeros(Nnewrois,2);
newroiVcorr = zeros(Nnewrois,1);

for ci = 1:Nnewrois
    tempiminds = sub2ind(size(imnewroi), news2p.stat{ci}.ypix, news2p.stat{ci}.xpix);
imnewroi(tempiminds) = ci;
newroictr(ci,:)=double(news2p.stat{ci}.med);
newroiVcorr(ci) = mean(news2pvcorr(tempiminds));
end

fprintf('news2p: on-ROI Vcorr mean value = %.4f\n', mean(news2pvcorr(imnewroi>0)))
fprintf('news2p: off-ROI Vcorr mean value = %.4f\n', mean(news2pvcorr(imnewroi==0)))

validnewXYcells = min(abs(newroictr(:,2)-cumsum([0 Lx'])),[],2)>cropedgethr & ...
    min(abs(newroictr(:,1)-[0 unique(Ly')]),[],2)>cropedgethr;

newroipairoverlap = false(Nnewrois,Nnewrois);
for ci = 1:Nnewrois
    tempiminds = sub2ind(size(imnewroi), news2p.stat{ci}.ypix, news2p.stat{ci}.xpix);
    col = unique(imnewroi(tempiminds));
    newroipairoverlap(ci,col(col>ci))=true;
end
% distance + correlation: for pairs with centers within 5 pixels and
% correlation larger than 0.2, discard the roi with the smaller dynamic range
newroipairumdist = sqrt( (yumperpix*(newroictr(:,1)-newroictr(:,1)')).^2 + (xumperpix*(newroictr(:,2)-newroictr(:,2)')).^2 );
newroipaircorr = corr(news2p.F');
[r,c] = find(newroipairumdist<=5 & newroipaircorr>0.2);
roi1s = r(r<c);
roi2s = c(r<c);
[r,c] = find(newroipairoverlap);
roi1s = [roi1s; r];
roi2s = [roi2s; c];
newroiinds2discard = zeros(size(roi1s));
for iroi = 1:numel(roi1s)
    if newroiVcorr(roi1s(iroi)) < newroiVcorr(roi2s(iroi))
        newroiinds2discard(iroi) = roi1s(iroi);
    else
        newroiinds2discard(iroi) = roi2s(iroi);
    end
end
newrois2keep = true(Nnewrois,1);
newrois2keep(newroiinds2discard) = false;
toc

newXYcoords = newroictr;
figure; hold all
plot(newroipairumdist(newroipairumdist<30), newroipaircorr(newroipairumdist<30), '.', 'Color', [0 0 0 0.1])
yl = ylim;
plot([5 5], yl, 'r-')
plot([0 30], [0.2 0.2], 'r-')
xlabel('news2p ROI pairwise distance (pixels)')
ylabel('news2p ROI pairwise correlation (raw F)')
title([mousedate ' all ROIs'], 'interpreter', 'none')

X2p = [4 8 16];
thresh2p = [0.012 0.015 0.02];
figure('Position', [0 0 2000 800])
subplot(2,2,1)
imshow(news2p.ops.Vcorr*10)
title([mousedate ' Vcorr'], 'interpreter', 'none')
for ii = 1:3
corriscell = newrois2keep & validnewXYcells & newroiVcorr>thresh2p(ii);
tempimroi = zeros(news2p.ops.Ly, news2p.ops.Lx);
tempimroi(ismember(imnewroi, find(corriscell))) = 1;
subplot(2,2,ii+1)
imshow(tempimroi)
title(sprintf('Vcorr threshold %.4f N=%d', thresh2p(ii), nnz(corriscell)))
end
Vcorrthresh = 0.02;
%%
addpath(genpath('d:\Users\USER\Documents\MATLAB\npy_matlab'))
newiscellnpy = readNPY([newpath 'iscell.npy']);

if ~isequal(size(newiscellnpy), size(news2p.iscell))
    error('check suite2p directory')
end
if ~isequal(newiscellnpy, news2p.iscell)
    warning('iscell in npy and Fall.mat are different -- perhaps I forgot to savemat after maunual curation?')
end

% for manual curation, use 0.02
newiscell = newrois2keep & validnewXYcells & newroiVcorr>Vcorrthresh;
newiscellnpy(:,1) = newiscell;
newiscellnpy(:,2) = newroiVcorr;
disp(nnz(newiscell))

figure; hold all
plot(newroipairumdist((newiscell'&newiscell) & newroipairumdist<30), newroipaircorr((newiscell'&newiscell) & newroipairumdist<30), '.', 'Color', [0 0 0 0.1])
yl = ylim;
plot([5 5], yl, 'r-')
plot([0 30], [0.2 0.2], 'r-')
xlabel('news2p ROI pairwise distance (pixels)')
ylabel('news2p ROI pairwise correlation (raw F)')
title([mousedate ' new iscell ROIs'], 'interpreter', 'none')

writeNPY(newiscellnpy, [newpath 'iscell.npy']);

load(strcat(newpath, 'Fall.mat'));
iscell = newiscellnpy;
% save(strcat(newpath, 'Fall.mat'), 'Vcorrthresh', 'F', 'iscell', 'redcell', 'stat', 'Fneu', 'ops', 'spks', '-v7.3')
%% online ROIs
% mean Vcorr per ROIs
Nonrois = numel(cat(2,online.stat));
onroictr = zeros(Nonrois,2);
onplane = zeros(Nonrois,1);

onlineVcorr = zeros(online(1).ops.Ly, online(1).ops.Lx, numel(online));
onlinemeanimg = zeros(online(1).ops.Ly, online(1).ops.Lx, numel(online));
onlinemaxproj = zeros(online(1).ops.Ly, online(1).ops.Lx, numel(online));

imonroi0 = zeros(online(1).ops.Ly, online(1).ops.Lx, numel(online));
onroiVcorr0 = zeros(Nonrois,1);
onroimeanimg0 = zeros(Nonrois,1);
onroimaxproj0 = zeros(Nonrois,1);
cnt = 0;
for iplane = 1:numel(online)
    % note, ops.Vcorr is not the same size as the input image!
    tempxinds = online(iplane).ops.xrange(1)+1:online(iplane).ops.xrange(2);
    tempyinds = online(iplane).ops.yrange(1)+1:online(iplane).ops.yrange(2);
    onlineVcorr(tempyinds,tempxinds,iplane) = online(iplane).ops.Vcorr;
    onlinemeanimg(:,:,iplane) = online(iplane).ops.meanImg;
    onlinemaxproj(tempyinds,tempxinds,iplane) = online(iplane).ops.max_proj;
    for ci = 1:numel(online(iplane).stat)
        cnt = cnt+1;
        onroictr(cnt,:)=double(online(iplane).stat{ci}.med);
        onplane(cnt) = iplane;
        tempiminds = sub2ind(size(imonroi0), online(iplane).stat{ci}.ypix+1, ...
            online(iplane).stat{ci}.xpix+1, iplane*ones(size(online(iplane).stat{ci}.xpix)));
        imonroi0(tempiminds) = ci;
        onroiVcorr0(cnt) = mean(onlineVcorr(tempiminds));
        onroimeanimg0(cnt) = mean(onlinemeanimg(tempiminds));
        onroimaxproj0(cnt) = mean(onlinemaxproj(tempiminds));
    end
end

fprintf('online: on-ROI Vcorr mean value = %.4f\n', mean(onlineVcorr(imonroi0>0)))
fprintf('online: on-ROI Vcorr mean value = %.4f\n', mean(onlineVcorr(imonroi0==0)))
%% 'iscell' criterion for online processed data
onroipairumdist = sqrt( (yumperpix*(onroictr(:,1)-onroictr(:,1)')).^2 + (xumperpix*(onroictr(:,2)-onroictr(:,2)')).^2);
onroipaircorr = corr(online.F');
[r,c] = find(onroipairumdist<=5 & onroipaircorr>0.2);
roi1s = r(r<c);
roi2s = c(r<c);
onroiinds2discard = zeros(size(roi1s));
for iroi = 1:numel(roi1s)
    if diff(prctile(online.F(roi1s(iroi),:), [1 99])) < ...
            diff(prctile(online.F(roi2s(iroi),:), [1 99]))
        onroiinds2discard(iroi) = roi1s(iroi);
    else
        onroiinds2discard(iroi) = roi2s(iroi);
    end
end
onrois2keep = true(Nonrois,1);
onrois2keep(onroiinds2discard) = false;
toc

figure; hold all
plot(onroipairumdist(onroipairumdist<30), onroipaircorr(onroipairumdist<30), '.', 'Color', [0 0 0 0.1])
yl = ylim;
plot([7 7], yl, 'r-')
plot([0 30], [0.2 0.2], 'r-')
xlabel('online ROI pairwise distance (pixels)')
ylabel('online ROI pairwise correlation (raw F)')
title([mousedate ' all ROIs'], 'interpreter', 'none')

newiscell = newrois2keep & validnewXYcells & newroiVcorr>Vcorrthresh;
oniscell = onrois2keep & onroiVcorr0>Vcorrthresh; % & validonXYcells
%% CHECK THAT THERE IS NO SYSTEMATIC SHIFT BETWEEN ONLINE AND OFFLINE ROI POSITIONS
figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', [mousedate ' online (r+) vs offline (bx) iscell ROIs'], 'edgecolor', 'none', 'interpreter', 'none')
hold all
temponrois = oniscell(isneuron);
plot(neuronXYcoords(temponrois,2), neuronXYcoords(temponrois,1), 'r+')
tempnewrois = newiscell;
plot(newXYcoords(tempnewrois,2), newXYcoords(tempnewrois,1), 'bx')
axis(double(news2p.ops.Lx)/2*[0 1 0 1])
box on
set(gca, 'YDir', 'reverse')
% make sure they're on the same plane
umdistthresh = 20;
disp('NOTE, DISTANCE IS OFF BY A FACTOR OF 2 IN EVERY SESSION')
disp('affected variables: onoffroiumdist, targetroiXYdist, offlineHoloRequest.target2celldist, targetroiumdist, minholodist')
onoffroiumdist = 2*sqrt( (yumperpix*(neuronXYcoords(:,1)-newXYcoords(:,1)')).^2 + ...
    (xumperpix*(neuronXYcoords(:,2)-newXYcoords(:,2)')).^2);
[mv,on2nearestoff]=nanmin(onoffroiumdist,[],2);

on2offinds = sub2ind(size(onoffroiumdist), (1:length(on2nearestoff))', on2nearestoff);
temp = neuronXYcoords(:,2)-newXYcoords(:,2)';
onoffroixshift = temp(on2offinds);
onoffroixshift(mv>=umdistthresh) = NaN;
temp = neuronXYcoords(:,1)-newXYcoords(:,1)';
onoffroiyshift = temp(on2offinds);
onoffroiyshift(mv>=umdistthresh) = NaN;

figure('Position', [100 100 1200 300])
annotation('textbox', [0.1 0.92 0.9 0.1], 'string', mousedate, 'edgecolor', 'none', 'interpreter', 'none')
subplot(1,3,1); 
hold all
histogram(mv, 'BinWidth', 1, 'Normalization', 'cdf')
yl = [0 1]; yl = ylim;
plot(umdistthresh*[1 1], yl, 'r-')
ylim(yl)
xlabel('distance to nearest offline ROI (pixels)')
title('per online ROI')
subplot(1,3,2); 
histogram(onoffroixshift, 'BinWidth', 1)
xlim([-10 10])
xlabel('X-shift (pixels)')
title(sprintf('mode %d', mode(onoffroixshift)))
subplot(1,3,3); 
histogram(onoffroiyshift, 'BinWidth', 1)
xlim([-10 10])
xlabel('Y-shift (pixels)')
title(sprintf('mode %d', mode(onoffroiyshift)))

%% offline_params: compare and align with presuite2p_params from check_before_s2pmeso
%% 온라인이랑 오프라인 상에서 tiff 파일 순서를 정확히 일치시켜야함. 
suite2pfilelist = news2p.ops.filelist;
suite2pfilelist(suite2pfilelist=='\')='/';
temp2p= cell(size(suite2pfilelist, 1), 1);
for f = 2:size(suite2pfilelist,1)%여기서부터. file00011 때문에 2부터 시작
 dotindex = strfind(suite2pfilelist(f,:), '.tif');
 temps2pfn = suite2pfilelist(f,1:dotindex+3);
 filepts = strsplit(temps2pfn, '/');
 temp2p{f}=fullfile(filepts{2:end-1});
end
temp2p = cellfun(@(x) char(x), temp2p, 'UniformOutput', false);
 dotindex = strfind(suite2pfilelist(1,:), '.tif');
 temps2pfn = suite2pfilelist(1,1:dotindex-1);
 filepts = strsplit(temps2pfn, '/');
 temp2p{1}=fullfile(filepts{2:end});
suite2pfilelist=unique(temp2p,'stable');
%여기까지는 로드하는 데이터가 tiff인 경우에만 해당됨
disp(suite2pfilelist)
frameind=[];

onlines2pparams = load([onlinepath 'presuite2p_params.mat']);
onlines2ph5folders = cellstr(expList(onlines2pparams.s2ph5folderind,:));

for f = 1:size(onlines2ph5folders,1)
whichfolder(f) = f;
fileindex(f) = f;
end

% numframess2pfile = zeros(size(onlines2ph5folders,1)*size(tiffs,1),1);
tic
% filenamesplitter ='h5file';
numframess2pfile = [];
for f = 1:size(suite2pfilelist,1)%여기에서 전체 tiff파일의 프레임 수를 로드하기,offline 경로
    import ScanImageTiffReader.ScanImageTiffReader.*;
    % dotindex = strfind(suite2pfilelist(f,:), '.h5');%이 부분은
    % fnsplit = strsplit(suite2pfilelist(f,1:dotindex-1), '/'); 
    % tiff인 경우 주석처리하고, h5인 경우 주석 풀고 밑 줄을 주석처리하면 됨
    fnsplit = strsplit(suite2pfilelist{f}, '\');
    tiffile=['D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\' fnsplit(end) '\'];
    tiffile=fullfile(tiffile{:});
    tiffile(tiffile=='/')='\';
    tiffs = dir(fullfile(tiffile, '*.tif'));

numframess2p = zeros(size(tiffs, 1), 1);
    if contains(tiffile, '\file')
    tiff_path = fullfile('D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106', 'file00011.tif');
    reader=ScanImageTiffReader(tiff_path);
    desc=reader.descriptions();
    numframess2p = size(desc,1);
    else
    for i = 1: size(tiffs,1)
    tiff_path = fullfile(tiffs(i).folder, tiffs(i).name);
    reader=ScanImageTiffReader(tiff_path);
    desc=reader.descriptions();
    numframess2p(i) = size(desc,1);
    end
    end
    numframess2pfile= [numframess2pfile; numframess2p];
    if any(strcmp(fnsplit{end}, onlines2ph5folders))
    startid(f) = (length(numframess2pfile)-length(numframess2p)+1);
    endid(f)=length(numframess2pfile);
    end
end
toc
startid= startid(startid ~= 0);
endid= endid(endid ~= 0);
% framindx=unique(framindx);
%% 
for f = 1:numel(allfoldernames)
    fileinds = fileindex(whichfolder==f);
    if ~isequal(fileinds, sort(fileinds))
        error('%s files not ordered correctly', allfoldernames{f})
    end
    if ~all(diff(fileinds)==1)
        warning('%s files non-consecutive indices', allfoldernames{f})
    end
end

fprintf('%d files with no vis index (home folder, spontaneous activity)\n', nnz(whichfolder==0))

tiffheader = imfinfo([refpath reftif]);
hSIh = tiffheader(1).Software;
hSIh = regexp(splitlines(hSIh), ' = ', 'split');
for n=1:length(hSIh)
    if strfind(hSIh{n}{1}, 'SI.hChannels.channelSave')
        nch = n;
        channelssaved = str2num(hSIh{n}{2});
    end
end
numchannels = numel(channelssaved);

numtimepointss2pfile = numframess2pfile/numchannels;
%% match online offline traces
% for each ONLINE ROI, find the OFFLINE ROI that is the most correlated.
% offline data should include the timepoints in the online data

off2ontimeinds = [];
for f = 1:numel(onlines2ph5folders)
starttimeind = [1; cumsum(numtimepointss2pfile(1:end-1))+1];
endtimeind = cumsum(numtimepointss2pfile);
off2ontimeinds = cat(2, off2ontimeinds, starttimeind(startid(f)):endtimeind(endid(f)));%여기에서 온라인 게 몇번째인지지 지정하는 것
end

if ~(length(off2ontimeinds)==size(online.F,2))
    error('error in matching time indices of online and offline processing')
end

onlineF = cat(1,online.F);
onoffcorrF = corr(onlineF', news2p.F(:,off2ontimeinds)');
if ~isequal(size(onoffroiumdist), size(onoffcorrF))
    error('check size mismatch')
end

[maxon2offcorrF1, imaxon2offcorrF] = max(onoffcorrF, [], 2);
[minon2offdist, iminon2offdist] = min(onoffroiumdist, [], 2);
fprintf('%.2f%% match in max F correlation and minimum distance online to offline\n', 100*mean(imaxon2offcorrF==iminon2offdist))

[maxoff2oncorrF, imaxoff2oncorrF] = max(onoffcorrF, [], 1);
[minoff2ondist, iminoff2ondist] = min(onoffroiumdist, [], 1);
fprintf('%.2f%% match in max F correlation and minimum distance offline to online\n', 100*mean(imaxoff2oncorrF==iminoff2ondist))

d1 = onoffroiumdist(sub2ind(size(onoffroiumdist), (1:length(imaxon2offcorrF))', imaxon2offcorrF));
figure; plot(d1, maxon2offcorrF1, '.')
xlabel('Distance (\mum)')
ylabel('Max On-Off F Correlation')
title('For Each Online ROI')
xlim([0 100])

% save(strcat(pathpp, 'online2offline.mat'), 'allfoldernames', 'onoffroiumdist', 'off2ontimeinds', 'onoffcorrF')
pair_count = 0;
threshold_corrF = 0.2; % maxon2offcorrF의 임계값
threshold_distance = 50;
for i = 1:numel(maxon2offcorrF1)
    if maxon2offcorrF1(i) > threshold_corrF && d(i) < threshold_distance
            pair_count = pair_count + 1;
    end
end
fprintf('maxon2offcorrF1 > %.2f and d < %.2f 인 pair의 개수: %d\n', threshold_corrF, threshold_distance, pair_count);
% 
figure;
hold on;
plot(d1, maxon2offcorrF1, '.', 'DisplayName', 'Original Data');
plot(d, maxon2offcorrF, '.', 'DisplayName', 'New Data');
xlabel('Distance (\mum)');
ylabel('Max On-Off F Correlation');
title('Comparison of Max On-Off F Correlation');
legend('Location', 'best');
xlim([0 500]);
hold off;