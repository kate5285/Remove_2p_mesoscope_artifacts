% check that my suite2p pipeline worked well by comparing with 
% online processing results (Mora's script) processed

% according to the Matlab suite2p github site, https://github.com/cortex-lab/Suite2P
% xpix, ypix: x and y indices of pixels belonging to this max. These index into the valid part of the image (defined by ops.yrange, ops.xrange).
% med: median of y and x pixels in the ROI (indices for the valid part of the image, defined by ops.yrange, ops.xrange).

clear all; %close all; clc
%% this first part is user-interactive: must be edited for each session
% takes about ~10 minutes. it's possible to skip this step if done in savio already
% mousedate = 'MU31_2/230111/';
mousedate = 'MU31_2/230106/';

disp(mousedate)

drivepath = '//shinlab/ShinLab/MesoHoloExpts/';
% mesoSIpath = 'D:\doyeon_kim\MesoHoloExpts_mesoholoexpts_scanimage_MU31_2\230106\h5file\';
mesoSIpath = [drivepath 'mesoholoexpts_scanimage/' mousedate];
onlinepath = [drivepath 'mesoholoexpts_scanimage/' mousedate 'ClosedLoop_justgreen/'];
vispath = [drivepath 'mesoholoexpts_scanimage/' mousedate 'visstiminfo/'];
path2p = [drivepath 'mesoholoexpts/' mousedate];
offlinepath = [drivepath 'mesoholoexpts/' mousedate 'suite2p/combined/'];
pathpp = [drivepath 'mesoholoexpts_postprocessed/' mousedate];

offline = load([offlinepath 'Fall.mat']);

cd([mesoSIpath 'visstiminfo'])
ls([mesoSIpath 'ClosedLoop_justgreen'])
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
%ls(mesoSIpath)
dircontents = dir(mesoSIpath);
mesoSIfolders = {dircontents([dircontents.isdir]).name};
disp(mesoSIfolders')
open mesoSIfolders

compute_numframes = true;

% input('USER NEEDS TO INPUT FOLDER NAMES FOR EVERY SESSION IN ORDER MATCHING EXPTIDN')
% %MU31_2/221227/
% foldername = {'101','retinotopy1','staticgratings','staticICtxi0','staticgratings12','staticICtxi1','SGholo'};
% %MU31_1/221227/ 
% foldername = {'101', 'retinotopy0', 'staticgratings', 'staticICtxi0', ...
%     'staticICtxi1', 'staticgratings12', 'stimtest', 'SGholo'};
% %MU31_2/230106/
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

%% vis files
exptids = cell(size(allfoldernames));
%nexpts = cell(size(foldername));
vis = struct();
for ii= 1:numel(allfoldernames)
    if strcmp(nexpts{ii}, 'x')
        if isletter(allfoldernames{ii}(1))
            exptids{ii} = allfoldernames{ii};
        else
            exptids{ii} = ['x' allfoldernames{ii}];
        end
    else
        vf = find(contains(Vordered, nexpts{ii}));
        load(strcat(vispath, Vordered{vf}))
        visvarname = who('-file',strcat(vispath, Vordered{vf}));
        visvar = eval(visvarname{1});
        
        exptids{ii} = visvar.exptid;
        if ~strcmp(nexpts{ii}, visvar.nexp)
            error('mismatch in nexp: check')
        end
        fprintf('%s matched to folder %s\n', Vordered{vf}, allfoldernames{ii})
        input('')
        
        exptidn = strcat(exptids{ii}, '_', nexpts{ii}); %exptidn은 뭐지
        %disp(exptidn)
        vis.(exptidn) = visvar;
        
        if ~isfield(vis.(exptidn), 'numtrials') || ~isfield(vis.(exptidn), 'durvisstim') || ~isfield(vis.(exptidn), 'iti')
            switch exptids{ii}
                case 'retinotopy'
                    vis.(exptidn).numtrials = size(vis.(exptidn).locinds, 1);
                    vis.(exptidn).durvisstim = numel(vis.(exptidn).orientations) * vis.(exptidn).nCycles/vis.(exptidn).tFreq;
                    vis.(exptidn).iti = vis.(exptidn).durvisstim * vis.(exptidn).ratio;
                case 'RFcircle'
                    vis.(exptidn).numtrials = size(vis.(exptidn).RFinds, 1);
                    vis.(exptidn).durvisstim = numel(vis.(exptidn).orientations) * vis.(exptidn).nCycles/vis.(exptidn).tFreq;
                    vis.(exptidn).iti = vis.(exptidn).durvisstim * vis.(exptidn).ratio;
                case {'gratings', 'staticgratings'}
                    if isfield(vis.(exptidn), 'durdrifting')
                        vis.(exptidn).durvisstim = vis.(exptidn).durdrifting;
                    else
                        vis.(exptidn).durvisstim = vis.(exptidn).durstatic;
                    end
                case 'blankICtxi'
                    continue
                otherwise
                    error('numtrials/durvisstim/iti must be specified')
            end
        end
        if strcmp(exptids{ii}, 'wheel')
            vis.(exptidn).durvisstim = vis.(exptidn).nCycles/vis.(exptidn).tFreq;
        end
        
        if ~isfield(vis.(exptidn), 'Telapsed') && isfield(vis.(exptidn), 'Tstartvisstim') && isfield(vis.(exptidn), 'Tendvisstim')
            vis.(exptidn).Telapsed = zeros(2*vis.(exptidn).numtrials, 1);
            vis.(exptidn).Telapsed(1:2:end-1) = vis.(exptidn).Tstartvisstim;
            vis.(exptidn).Telapsed(2:2:end) = vis.(exptidn).Tendvisstim;
        end
        
    end
end

% save(strcat(pathpp, 'vis_params.mat'), 'exptids', 'nexpts', 'vis')

%%
% RESTRICT FOV IN OFFLINE SUITE2P?'
restrictofffov = false; % if false, offxlb, offxub, offylb and offyub will be ignored
offxlb = 450; offxub = 2450;
offylb = 0; offyub = size(offline.ops.meanImg,1);

[holoreqfile,holoreqpath] = uigetfile([onlinepath '*.mat'], 'holoRequest_MU31_2_230106_staticICtxi_001');%이게 맞는지 확인하기
load([holoreqpath holoreqfile])
fullnpix = [sum(nxpix),mean(nypix)];
[fullxpix,fullypix] = meshgrid(linspace(1,fullnpix(1),fullnpix(1)),...
    linspace(1,fullnpix(2),fullnpix(2)));
fullxsize = sum(xsize); %um
fullysize = mean(ysize); %um
fullxcenter = (xcenter(1)-xsize(1)/2) + ...
    abs((xcenter(1)-xsize(1)/2) - (xcenter(end)+xsize(end)/2))/2; %um
fullycenter = mean(ycenter); %um
xumperpix = fullxsize/size(offline.ops.meanImg,2);
yumperpix = fullysize/size(offline.ops.meanImg,1);

%% EXCLUDE NEURONS AT THE EDGE OF EACH STRIPE
cropedgethr = 0; % 0 to not restrict any cells %여기 바꿈

% copy pasted from mesoscope_json_from_scanimage
[reftif, refpath] = uigetfile('\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts_scanimage\MU31_2\230106\sizecircleC\file_00210.tif');
% 'M:\mesoholoexpts_scanimage\HS_Chrome2f_3\220923\ref\file_00001.tif';
fname = [refpath reftif];

header = imfinfo(fname);
% stack = loadFramesBuff(fname,1,1,1);
% % get relevant infomation for TIFF header
artist_info     = header(1).Artist;
% % retrieve ScanImage ROIs information from json-encoded string 
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

if ~(sum(Lx)==offline.ops.Lx && isequal(unique(Ly), offline.ops.Ly))
    error('unexpected Lx/Ly error')
end

%%
load([onlinepath 'online_params.mat'])
% load(sprintf('%sholoRequest_%s_%s.mat', onlinepath, animalid, dstr))
online = load([onlinepath '/suite2p/plane0/Fall.mat']);

offstat = cell2mat(offline.stat);
offmed = double( cat(1, offstat.med) );
% offplane = double(1,1,1,1,1,1);%여기를 바꿈

% 1 + cat(1, offstat.iplane)
%%%%%%%%%%%%%%%% MEAN IMAGE %%%%%%%%%%%%%%%%
fs = 16;
figure
imshow(online.ops.meanImg/200)
%caxis([0 100])
%colorbar
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', 'online suite2p mean img', 'EdgeColor', 'None', 'FontSize', fs)

figure
imshow(offline.ops.meanImg/200)
% caxis([0 300])
% title('offline suite2p mean img')
colormap gray
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', [mousedate ' offline suite2p mean img'], 'EdgeColor', 'none', 'FontSize', fs, 'interpreter', 'none')

% figure; imagesc(offline.ops.meanImg_chan2)
% caxis([0 500])


%%%%%%%%%%%%%%%% MEAN IMAGE CH2 %%%%%%%%%%%%%%%%
if offline.ops.nchannels>1
figure; imshow(offline.ops.meanImg_chan2/100)
% caxis([0 300])
% title('offline suite2p mean img')
colormap gray
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', [mousedate ' offline suite2p red channel mean img'], 'EdgeColor', 'None', 'FontSize', fs, 'interpreter', 'none')

% figure; imagesc(offline.ops.meanImg_chan2)
% caxis([0 500])
end
%}

%%%%%%%%%%%%%%%% Vcorr IMAGE %%%%%%%%%%%%%%%%
fs = 16;
% figure; imshow(offline.ops.Vcorr*4)
% % caxis([0 300])
% % title('offline suite2p mean img')
% colormap gray
% annotation('textbox', [0.1 0.9 0.9 0.1], 'String', 'offline suite2p Vcorr', 'EdgeColor', 'None', 'FontSize', fs)

figure
imshow(offline.ops.Vcorr(y_start:y_end,x_start:x_end)*4)
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', [mousedate ' offline suite2p Vcorr: for direct comparison'], 'EdgeColor', 'None', 'FontSize', fs, 'interpreter', 'none')

figure
imshow(online.ops.Vcorr*4)
annotation('textbox', [0.1 0.9 0.9 0.1], 'String', [mousedate ' online suite2p Vcorr'], 'EdgeColor', 'None', 'FontSize', fs, 'interpreter', 'none')

% %% 
% figure; hold all
% plot(offline.ops.xoff)
% plot(offline.ops.yoff)

%% 'iscell' criterion for pre-manual curation
tic
% mean Vcorr per ROIs
Noffrois = numel(offline.stat);
imoffroi = zeros(size(offline.ops.Vcorr));%수정
offroictr = zeros(Noffrois,2);
offroiVcorr = zeros(Noffrois,1);

for ci = 1:Noffrois
    valid_indicesx = find(offline.stat{ci}.xpix <= 2168);%수정
    valid_indicesy = find(offline.stat{ci}.ypix <= 1046);
    
    % Check if valid_indicesx and valid_indicesy have the same length
    if length(valid_indicesx) ~= length(valid_indicesy)
        continue; % Skip this iteration if lengths don't match
    end
    
    offline.stat{ci}.ypix = offline.stat{ci}.ypix(valid_indicesy);
    offline.stat{ci}.xpix = offline.stat{ci}.xpix(valid_indicesx);

    tempiminds = sub2ind(size(offline.ops.Vcorr), offline.stat{ci}.ypix, offline.stat{ci}.xpix);
    imoffroi(tempiminds) = ci;
    offroictr(ci, :) = double(offline.stat{ci}.med);
    offroiVcorr(ci) = mean(offline.ops.Vcorr(tempiminds));
end
fprintf('offline: on-ROI Vcorr mean value = %.4f\n', mean(offline.ops.Vcorr(imoffroi>0)))
fprintf('offline: off-ROI Vcorr mean value = %.4f\n', mean(offline.ops.Vcorr(imoffroi==0)))
% disp(median(offline.ops.Vcorr(imoffroi==0)))

% validXYcells: the point is to get rid of stim artifact.
% one option is to use the same criterion as online prcessing.
validoffXYcells = min(abs(offroictr(:,2)-cumsum([0 Lx'])),[],2)>cropedgethr & ...
    min(abs(offroictr(:,1)-[0 unique(Ly')]),[],2)>cropedgethr;
if restrictofffov
    validoffXYcells = validoffXYcells & offroictr(:,2)>offxlb & offroictr(:,2)<=offxub & offroictr(:,1)>offylb & offroictr(:,1)<=offyub;
end

% extra criteria: no overlapping pixels between ROIs
% pick out cells that have overlapping rois: fill up upper right triangle
offroipairoverlap = false(Noffrois,Noffrois);
for ci = 1:Noffrois
   valid_indicesx = find(offline.stat{ci}.xpix <= 2168);%수정
    valid_indicesy = find(offline.stat{ci}.ypix <= 1046);
    
    % Check if valid_indicesx and valid_indicesy have the same length
    if length(valid_indicesx) ~= length(valid_indicesy)
        continue; % Skip this iteration if lengths don't match
    end
    
    offline.stat{ci}.ypix = offline.stat{ci}.ypix(valid_indicesy);
    offline.stat{ci}.xpix = offline.stat{ci}.xpix(valid_indicesx);

    tempiminds = sub2ind(size(imoffroi), offline.stat{ci}.ypix, offline.stat{ci}.xpix);
    col = unique(imoffroi(tempiminds));
    offroipairoverlap(ci,col(col>ci))=true;
end

% distance + correlation: for pairs with centers within 5 pixels and
% correlation larger than 0.2, discard the roi with the smaller dynamic range
offroipairumdist = sqrt( (yumperpix*(offroictr(:,1)-offroictr(:,1)')).^2 + (xumperpix*(offroictr(:,2)-offroictr(:,2)')).^2 );
offroipaircorr = corr(offline.F');
[r,c] = find(offroipairumdist<=5 & offroipaircorr>0.2);
roi1s = r(r<c);
roi2s = c(r<c);
[r,c] = find(offroipairoverlap);
roi1s = [roi1s; r];
roi2s = [roi2s; c];
offroiinds2discard = zeros(size(roi1s));
for iroi = 1:numel(roi1s)
    if offroiVcorr(roi1s(iroi)) < offroiVcorr(roi2s(iroi))
        offroiinds2discard(iroi) = roi1s(iroi);
    else
        offroiinds2discard(iroi) = roi2s(iroi);
    end
end
offrois2keep = true(Noffrois,1);
offrois2keep(offroiinds2discard) = false;
toc

offXYcoords = offroictr;
% figure; histogram(offroiVcorr)
% figure; histogram(offroipaircorr(offroipairumdist>0 & offroipairumdist<5), 'BinWidth', 0.05)

figure; hold all
plot(offroipairumdist(offroipairumdist<30), offroipaircorr(offroipairumdist<30), '.', 'Color', [0 0 0 0.1])
yl = ylim;
plot([5 5], yl, 'r-')
plot([0 30], [0.2 0.2], 'r-')
xlabel('offline ROI pairwise distance (pixels)')
ylabel('offline ROI pairwise correlation (raw F)')
title([mousedate ' all ROIs'], 'interpreter', 'none')

X2p = [4 8 16];
thresh2p = [0.012 0.015 0.02];
figure('Position', [0 0 2000 800])
subplot(2,2,1)
imshow(offline.ops.Vcorr*10)
title([mousedate ' Vcorr'], 'interpreter', 'none')
for ii = 1:3
% corriscell = offroiVcorr>X2p(ii)*mean(offline.ops.Vcorr(imoffroi==0));
corriscell = offrois2keep & validoffXYcells & offroiVcorr>thresh2p(ii);
tempimroi = zeros(offline.ops.Ly, offline.ops.Lx);
tempimroi(ismember(imoffroi, find(corriscell))) = 1;
subplot(2,2,ii+1)
imshow(tempimroi)
title(sprintf('Vcorr threshold %.4f N=%d', thresh2p(ii), nnz(corriscell)))
end

% figure; histogram(offroiVcorr(offrois2keep & validoffXYcells))
% xlabel('offroiVcorr')

%% USER INPUT Vcorrthresh
Vcorrthresh = 0.02; %input('USER INPUT Vcorrthresh')
%% rewrite iscell.npy https://github.com/kwikteam/npy-matlab
addpath(genpath('d:\Users\USER\Documents\MATLAB\npy_matlab'))
% a = rand(5,4,3);
% writeNPY(a, 'a.npy');
% b = readNPY('a.npy');
offiscellnpy = readNPY([offlinepath 'iscell.npy']);

if ~isequal(size(offiscellnpy), size(offline.iscell))
    error('check suite2p directory')
end
if ~isequal(offiscellnpy, offline.iscell)
    warning('iscell in npy and Fall.mat are different -- perhaps I forgot to savemat after maunual curation?')
end

% for manual curation, use 0.02
offiscell = offrois2keep & validoffXYcells & offroiVcorr>Vcorrthresh;
offiscellnpy(:,1) = offiscell;
offiscellnpy(:,2) = offroiVcorr;
disp(nnz(offiscell))

figure; hold all
plot(offroipairumdist((offiscell'&offiscell) & offroipairumdist<30), offroipaircorr((offiscell'&offiscell) & offroipairumdist<30), '.', 'Color', [0 0 0 0.1])
yl = ylim;
plot([5 5], yl, 'r-')
plot([0 30], [0.2 0.2], 'r-')
xlabel('offline ROI pairwise distance (pixels)')
ylabel('offline ROI pairwise correlation (raw F)')
title([mousedate ' offiscell ROIs'], 'interpreter', 'none')


writeNPY(offiscellnpy, [offlinepath 'iscell.npy']);

load(strcat(offlinepath, 'Fall.mat')); % F, Fneu, iscell, ops, redcell, spks, stat
iscell = offiscellnpy;
save(strcat(offlinepath, 'Fall.mat'), 'Vcorrthresh', 'F', 'iscell', 'redcell', 'stat', 'Fneu', 'ops', 'spks', '-v7.3')


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

% validXYcells: the point is to get rid of stim artifact.
% one option is to use the same criterion as online prcessing.
% the most permissive is to start at x=11
% bound = 50;
% validonXYcells = onroictr(:,2)>xstart & onroictr(:,2)<=xend & ...
%     onroictr(:,1)>bound & onroictr(:,1)<=online(1).ops.Ly-bound;

% distance + correlation: for pairs with centers within 5 pixels and
% correlation larger than 0.2, discard the roi with the smaller dynamic range
onroipairumdist = sqrt( (yumperpix*(onroictr(:,1)-onroictr(:,1)')).^2 + (xumperpix*(onroictr(:,2)-onroictr(:,2)')).^2);
% onplanediff = onplane - onplane';
% onroipairumdist(onplanediff~=0)=NaN;

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

offiscell = offrois2keep & validoffXYcells & offroiVcorr>Vcorrthresh;
oniscell = onrois2keep & onroiVcorr0>Vcorrthresh; % & validonXYcells
% figure; histogram(onroiVcorr0); xlabel('online ROI Vcorr')

%% CHECK THAT THERE IS NO SYSTEMATIC SHIFT BETWEEN ONLINE AND OFFLINE ROI POSITIONS
figure('Position', [0 0 1000 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', [mousedate ' online (r+) vs offline (bx) iscell ROIs'], 'edgecolor', 'none', 'interpreter', 'none')
hold all
temponrois = oniscell(isneuron);
plot(neuronXYcoords(temponrois,2), neuronXYcoords(temponrois,1), 'r+')
tempoffrois = offiscell;
plot(offXYcoords(tempoffrois,2), offXYcoords(tempoffrois,1), 'bx')
axis(double(offline.ops.Lx)/2*[0 1 0 1])
box on
set(gca, 'YDir', 'reverse')
%     axis off

% for each ONLINE ROI, find the OFFLINE ROI that is nearest.
% make sure they're on the same plane
umdistthresh = 20;
disp('NOTE, DISTANCE IS OFF BY A FACTOR OF 2 IN EVERY SESSION')
disp('affected variables: onoffroiumdist, targetroiXYdist, offlineHoloRequest.target2celldist, targetroiumdist, minholodist')
onoffroiumdist = 2*sqrt( (yumperpix*(neuronXYcoords(:,1)-offXYcoords(:,1)')).^2 + ...
    (xumperpix*(neuronXYcoords(:,2)-offXYcoords(:,2)')).^2);
% onplanediff = onplane(isneuron) - offplane';
% onoffroiumdist(onplanediff~=0)=NaN;
[mv,on2nearestoff]=nanmin(onoffroiumdist,[],2);

%{
Zfactor = 3.5*2; % this should correspond to the axial-PSF/lateral-PSF * micrometer/pixel
fprintf('CHECK THAT axial-PSF/lateral-PSF RATIO X um/pixel IS %.1f\n', Zfactor)
offZcoords = planeZcoords(offplane)';
onoffXYZpixdist = sqrt( (neuronXYcoords(:,1)-offXYcoords(:,1)').^2 + ...
    (neuronXYcoords(:,2)-offXYcoords(:,2)').^2 + ...
    ((neuronXYcoords(:,3)-offZcoords')/Zfactor).^2 );
% figure; plot(onoffroipixdist, onoffroiXYZpixdist, '.')
% xlabel('On-Off XY dist')
% ylabel('On-Off XYZ dist')
% title('expecting a perfect diagonal')
onoffZdist = neuronXYcoords(:,3)-offZcoords';
if ~all(isnan(onoffroipixdist(onoffZdist~=0)))
    error('expecting onoffroipixdist to be NaN for ROIs in different planes')
end
%}

% on2offinds = sub2ind(size(onoffroiumdist), find(mv<umdistthresh), on2nearestoff(mv<umdistthresh));
on2offinds = sub2ind(size(onoffroiumdist), (1:length(on2nearestoff))', on2nearestoff);
temp = neuronXYcoords(:,2)-offXYcoords(:,2)';
onoffroixshift = temp(on2offinds);
onoffroixshift(mv>=umdistthresh) = NaN;
temp = neuronXYcoords(:,1)-offXYcoords(:,1)';
onoffroiyshift = temp(on2offinds);
onoffroiyshift(mv>=umdistthresh) = NaN;

% figure%('Position', [100 100 1200 300])
% subplot(2,2,1); 
% hold all
% histogram(mv, 'BinWidth', 1, 'Normalization', 'cdf')
% yl = [0 1]; yl = ylim;
% plot(umdistthresh*[1 1], yl, 'r-')
% ylim(yl)
% xlabel('distance to nearest offline ROI (per online ROI, pixels)')
% subplot(2,2,3); 
% histogram(onoffroixshift, 'BinWidth', 1)
% xlabel('X-shift (pixels)')
% subplot(2,2,4); 
% histogram(onoffroiyshift, 'BinWidth', 1)
% xlabel('Y-shift (pixels)')


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

% figure;
% % imshow(offline.ops.meanImg/300)
% imshow(10*offline.ops.Vcorr)
% hold on
% plot(offroictr(offiscell,2), offroictr(offiscell,1), 'r.')

%% offline_params: compare and align with presuite2p_params from check_before_s2pmeso
% 605번 줄에서 막힘, h5이라 이 부분 스킵할 것.
suite2pfilelist = offline.ops.filelist;
% if ~(numel(pres2p.tiffns)==size(suite2pfilelist,1))
%     error('suite2pfilelist should match presuite2p_params files')
% end
% tiffnind = zeros(size(suite2pfilelist,1),1);
whichfolder = zeros(size(suite2pfilelist,1),1);
fileindex = zeros(size(suite2pfilelist,1),1);
for f = 1:size(suite2pfilelist,1)
    dotindex = strfind(suite2pfilelist(f,:), '.h5');%h5이라 이 부분 스킵할 것.
    temps2pfn = suite2pfilelist(f,1:dotindex-1);
%     slashindex = strfind(temps2pfn, '/');
%     matchingpres2pfile = strcmp(tiffnfoldersplit(:,end), temps2pfn(slashindex(end-1)+1:slashindex(end)-1)) ...
%         & strcmp({pres2p.tiffns.name}', temps2pfn(slashindex(end)+1:end));
%     if nnz(matchingpres2pfile) ~= 1
%         error('tiffnind check failed')
%     end
%     tiffnind(f) = find(matchingpres2pfile);
%     % disp(suite2pfilelist(f,:))
%     % disp([pres2p.tiffns(tiffnind(f)).folder '/' pres2p.tiffns(tiffnind(f)).name])
    % fnsplit = strsplit(suite2pfilelist(f,1:dotindex-1), '_');%여기 주석처리함
    % fileindex(f) = str2num(fnsplit{end});
    temps2pfn(temps2pfn=='\') = '/';
    fileparts = strsplit(temps2pfn, '/');%수정
    folderind = find(strcmp(allfoldernames, fileparts{end}));%수정
    if numel(folderind)>1
        error('more than one corresponding vis file? check')
    end
    if isempty(folderind)
        whichfolder(f) = 0; % blank file with no vis index
    else
        whichfolder(f) = folderind;
    end    
end
%% 이 부분은 permission denied됨. (synology에 write file할 권리,666번 줄)
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

% takes about ~10 minutes. it's possible to skip this step if done in savio already
if exist(strcat(pathpp, 'offline_params.mat'))
    load(strcat(pathpp, 'offline_params.mat'), 'numframess2pfile')
else
tic
numframess2pfile = zeros(size(suite2pfilelist,1),1);
filenamesplitter = mousedate;
for f = 1:size(suite2pfilelist,1)
    import ScanImageTiffReader.ScanImageTiffReader.*;
    fileparts = strsplit(suite2pfilelist(f,:), filenamesplitter);
    dotindex = strfind(fileparts{end}, '.tif');
    filename = fileparts{end}(1:dotindex+3);
    tiffile = [mesoSIpath filename];
    
    reader=ScanImageTiffReader(tiffile);
    desc=reader.descriptions();
    numframess2pfile(f) = size(desc,1);
end
toc
end

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
Nplanes = double(offline.ops.nplanes);

numtimepointss2pfile = numframess2pfile/numchannels;
save(strcat(pathpp, 'offline_params.mat'), 'allfoldernames', 'whichfolder', 'exptids', 'nexpts', 'vis', ...
    'Nplanes', 'numchannels', 'numtimepointss2pfile', 'numframess2pfile', 'holoreqfile', 'holoreqpath')

%'tiffnind', 
%% match online offline traces
% for each ONLINE ROI, find the OFFLINE ROI that is the most correlated.
% offline data should include the timepoints in the online data
% off2ontimeinds 
onlines2pparams = load([onlinepath 'presuite2p_params.mat']);
onlines2ph5folders = cellstr(expList(onlines2pparams.s2ph5folderind,:));
off2ontimeinds = [];
for f = 1:numel(onlines2ph5folders)
% whichfolder = input('index of online-processed folder\n');
folderind = find(strcmp(onlines2ph5folders{f}, allfoldernames));
startfileind = find(whichfolder==folderind, 1, 'first');
endfileind = find(whichfolder==folderind, 1, 'last');
starttimeind = [1; cumsum(numtimepointss2pfile(1:end-1))+1];
endtimeind = cumsum(numtimepointss2pfile);

off2ontimeinds = cat(2, off2ontimeinds, starttimeind(startfileind):endtimeind(endfileind));%startfileind와 endfileind도 안 해도 될 것 같음.
end

if ~( length(off2ontimeinds)==size(online.F,2) ) %온라인 파일 상으로는 h5로 전환한 것이 3개.closed loop green안의 F와 안 맞는 것 같음
    %
    error('error in matching time indices of online and offline processing')
end

onlineF = cat(1,online.F);
onoffcorrF = corr(onlineF', offline.F(:,off2ontimeinds)');
if ~isequal(size(onoffroiumdist), size(onoffcorrF))
    error('check size mismatch')
end

[maxon2offcorrF, imaxon2offcorrF] = max(onoffcorrF, [], 2);
[minon2offdist, iminon2offdist] = min(onoffroiumdist, [], 2);
fprintf('%.2f%% match in max F correlation and minimum distance online to offline\n', 100*mean(imaxon2offcorrF==iminon2offdist))

[maxoff2oncorrF, imaxoff2oncorrF] = max(onoffcorrF, [], 1);
[minoff2ondist, iminoff2ondist] = min(onoffroiumdist, [], 1);
fprintf('%.2f%% match in max F correlation and minimum distance offline to online\n', 100*mean(imaxoff2oncorrF==iminoff2ondist))

d = onoffroiumdist(sub2ind(size(onoffroiumdist), (1:length(imaxon2offcorrF))', imaxon2offcorrF));
figure; plot(d, maxon2offcorrF, '.')
xlabel('Distance (\mum)')
ylabel('Max On-Off F Correlation')
title('For Each Online ROI')
xlim([0 100])

save(strcat(pathpp, 'online2offline.mat'), 'allfoldernames', 'onoffroiumdist', 'off2ontimeinds', 'onoffcorrF')
