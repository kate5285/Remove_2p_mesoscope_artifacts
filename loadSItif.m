% generate a avg image from a raw SI tif file
% check that this matches with the suite2p meanimg, i.e., that it is not flipped in any dimension

% How can I select multiple points using the Data Cursor and export the coordinates to the MATLAB workspace?
% https://www.mathworks.com/matlabcentral/answers/94353-how-can-i-select-multiple-points-using-the-data-cursor-and-export-the-coordinates-to-the-matlab-work
% 1. Generate figure.
% 2. Click the Data Cursor button on the toolbar of the generated figure.
% 3. Click any point of your choice on the line in the figure.
% 4. While pressing the Alt key, repeat step 3 as many times as you like until you have selected your desired set of points.
% 5. Right-click (or control-click if you are on a Mac) anywhere on the figure, and select the 'Export Cursor Data to Workspace...' option from the context menu.
% 6. Accept the default variable name, "cursor_info", and click "OK".
% 7. Type "cursor_info.Position" at the MATLAB command prompt and hit "Enter".

function [meanimg, maxprojimg]=loadSItif(tiffn, pltopt)
    addpath(genpath('d:\Users\USER\Documents\MATLAB'))
% tiffn='D:\doyeon_kim\MesoHoloExpts_scanimage_MU31_2_ICholo\230106\file_00001.tif'    

header = imfinfo(tiffn);
hSIh = header(1).Software;
hSIh = regexp(splitlines(hSIh), ' = ', 'split');
for n=1:length(hSIh)
	if strfind(hSIh{n}{1}, 'SI.hRoiManager.scanVolumeRate')
        fprintf('%d %s: %s\n', n, hSIh{n}{1}, hSIh{n}{2})
		fs = str2num(hSIh{n}{2});
	end
    if strfind(hSIh{n}{1}, 'SI.hFastZ.userZs')
        fprintf('%d %s: %s\n', n, hSIh{n}{1}, hSIh{n}{2})
		zs = str2num(hSIh{n}{2});
		nplanes = numel(zs);
    end
    if strfind(hSIh{n}{1}, 'SI.hChannels.channelSave')
        fprintf('%d %s: %s\n', n, hSIh{n}{1}, hSIh{n}{2})
        nchannels = numel(str2num(hSIh{n}{2}));
    end
	
end
artist_info     = header(1).Artist;
artist_info = artist_info(1:find(artist_info == '}', 1, 'last'));
artist = jsondecode(artist_info);

%%
si_rois = artist.RoiGroups.imagingRoiGroup.rois;
% get ROIs dnsions for each z-plane
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

% deduce flyback frames from most filled z-plane.아래는 총 프레임 수 구하는 함수.
function nFrames = nFramesTiff(tiffn)
tiff = Tiff(tiffn, 'r');
    closeTiff = onCleanup(@() close(tiff));
    
    nFrames = 0;
    while ~lastDirectory(tiff)
        nFrames = nFrames + 1;
        nextDirectory(tiff);
    end
end
%for selecting few frames, use firstIdx and lastIdx like the ones below.
stack = loadFramesBuff(tiffn,1,nFramesTiff(tiffn));

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
import ScanImageTiffReader.ScanImageTiffReader.*
%sipath = 'D:\doyeon_kim\ScanImageTiffReader';
%tiffns(1).name = 'file_00001.tif';

%f0=ScanImageTiffReader([sipath tiffns(1).name]);
%f0data = f0.data();
%f0mean = mean(f0data,3);

imfile = ScanImageTiffReader();
imfile.open(tiffn);
imdata = imfile.data();

immean = cell(nchannels, nplanes);
immaxproj = cell(nchannels, nplanes);
for ich = 1:nchannels
    for iplane = 1:nplanes
        zinds = (ich-1)+(iplane-1)*nchannels+1:nchannels*nplanes:size(imdata,3);
immean{ich,iplane} = squeeze(mean(imdata(:,:,zinds),3));
immaxproj{ich,iplane} = squeeze(max(imdata(:,:,zinds),[],3));
    end
end

meanimg = cell(nchannels, nplanes);
maxprojimg = cell(nchannels, nplanes);
for ich = 1:nchannels
    for iplane = 1:nplanes
        meanimg{ich, iplane} = zeros(sum(Ly), max(Lx));
        maxprojimg{ich, iplane} = zeros(sum(Ly), max(Lx));
        for istrip=1:numel(Ly)
            disp(numel(Ly))
            meanimg{ich,iplane}(sum(Ly(1:istrip-1))+1:sum(Ly(1:istrip)),1:Lx(istrip))=immean{ich, iplane}(:, data.lines{istrip}+1)';
            maxprojimg{ich,iplane}(sum(Ly(1:istrip-1))+1:sum(Ly(1:istrip)),1:Lx(istrip))=immean{ich, iplane}(:, data.lines{istrip}+1)';
        end
    end
end


for ich = 1:nchannels
    for iplane = 1:nplanes
figure; hold all
imagesc(meanimg{ich,iplane}); caxis([-20 50]); colorbar
colormap gray
set(gca, 'YDir', 'reverse')
title(sprintf('Plane%d PMT%d', ich-1, iplane-1))
    end
end
% plot(xynew(:,1), xynew(:,2), 'rx')


%%
%red_frames = stack(:,:,2:2:end);
%redness = zeros(1, size(red_frames, 3));
%for i = 1:size(red_frames, 3)
    %frame = red_frames(:,:,i);
    %redness(i) = mean2(frame);
%end

% x축에 프레임 번호, y축에 redness 출력
    %figure;
    %plot(1:numel(redness), redness);
    %xlabel('프레임 번호');
    %ylabel('redness');
    %title('프레임별 Redness');
meanimg_plane2 = meanimg{2, 1};
mean_vertical = mean(meanimg_plane2, 1);
x = 1:size(mean_vertical, 2);
figure;
plot(x, mean_vertical);
xlabel('location along plane width');
ylabel('mean redness');
title('(red channel) mean redness along plane x axis');
end
