%% run after running check_onlineVSoffline if necessary
% for MBOT45_007/220509, had a unique problem where online to offline match
% was reasonble for all planes other than the third plane.
% for HS_EmxChrmine_3/220710, and MBOT55/220710 all planes shifted the same amount. 

load([onlinepath 'online_params.mat'])%, 'neuronXYZcoords')
neuronXYcoords0 = neuronXYcoords;

onoffroiumdist = 2*sqrt( (yumperpix*(neuronXYcoords(:,1)-newXYcoords(:,1)')).^2 + ...
    (xumperpix*(neuronXYcoords(:,2)-newXYcoords(:,2)')).^2);

d = onoffroiumdist(sub2ind(size(onoffroiumdist), (1:nnz(isneuron))', imaxon2offcorrF));
figure; plot(d, maxon2offcorrF, '.')
xlim([0 50])
xlabel('Distance (Pixels)')
ylabel('Max On-Off F Correlation')
title('For Each Online ROI')

% figure; histogram(d(maxon2offcorrF>0.5))

nonmatch = maxon2offcorrF>0.5;% & d>20;

figure; hold all
plot(neuronXYcoords(nonmatch,2), neuronXYcoords(nonmatch,1), 'r+')
plot(newXYcoords(imaxon2offcorrF(nonmatch),2), newXYcoords(imaxon2offcorrF(nonmatch),1), 'bx')

xydiff = newXYcoords(imaxon2offcorrF(nonmatch),:) - neuronXYcoords(nonmatch,:);

limaxes = true;
xl = [-10 10];
figure('Position', [100 100 1200 300])
subplot(1,3,1)
histogram(xydiff(:,1), 'BinWidth', 1)
title(sprintf('X off-on mode %d', mode(xydiff(:,1))))
if limaxes
xlim(xl)
end
subplot(1,3,2)
histogram(xydiff(:,2), 'BinWidth', 1)
if limaxes
xlim(xl)
end
title(sprintf('Y off-on mode %d', mode(xydiff(:,2))))
subplot(1,3,3)
histogram2(xydiff(:,1), xydiff(:,2), 'BinWidth', 1, 'DisplayStyle', 'tile')
if limaxes
axis([xl xl])
end
xlabel('X-shift')
ylabel('Y-shift')

% %% correct for the shift on the third plane
% offXYZcoords(offplane==3,:) = offXYZcoords(offplane==3,:) - mode(xyzdiff,1); % for MBOT45_007/220509
newXYcoords = newXYcoords - mode(xydiff,1); % for HS_EmxChrmine_3/220710
% clearvars offXYcoords
% neuronXYZcoords(onplane(isneuron)==3,:) = neuronXYZcoords(onplane(isneuron)==3,:) + mode(xyzdiff,1);
% FOLLOWING SECTIONS ARE COPY-PASTED FROM check_onlineVSoffline AND SHOULD BE IDENTICAL

%% CHECK THAT THERE IS NO SYSTEMATIC SHIFT BETWEEN ONLINE AND OFFLINE ROI POSITIONS
% onplane and offplane denotes which stripe in this code
figure('Position', [0 0 1800 1000])
annotation('textbox', [0.1 0.9 0.9 0.1], 'string', 'online (r+) vs offline (bx) iscell ROIs', 'edgecolor', 'none')
hold all
temponrois = oniscell(isneuron);
plot(neuronXYcoords(temponrois,2), neuronXYcoords(temponrois,1), 'r+')
tempoffrois = newiscell;
plot(newXYcoords(tempoffrois,2), newXYcoords(tempoffrois,1), 'bx')

% for each ONLINE ROI, find the OFFLINE ROI that is nearest.
% make sure they're on the same plane
umdistthresh = 10;
onoffroiumdist = 2*sqrt( (yumperpix*(neuronXYcoords(:,1)-newXYcoords(:,1)')).^2 + ...
    (xumperpix*(neuronXYcoords(:,2)-newXYcoords(:,2)')).^2);
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
temp = neuronXYcoords(:,2)-newXYcoords(:,2)';
onoffroixshift = temp(on2offinds);
onoffroixshift(mv>=umdistthresh) = NaN;
temp = neuronXYcoords(:,1)-newXYcoords(:,1)';
onoffroiyshift = temp(on2offinds);
onoffroiyshift(mv>=umdistthresh) = NaN;

figure%('Position', [100 100 1200 300])
subplot(2,2,1); 
hold all
histogram(mv, 'BinWidth', 1, 'Normalization', 'cdf')
yl = [0 1]; yl = ylim;
plot(umdistthresh*[1 1], yl, 'r-')
ylim(yl)
xlabel('distance to nearest offline ROI (per online ROI, pixels)')
subplot(2,2,3); 
histogram(onoffroixshift, 'BinWidth', 1)
xlabel('X-shift (pixels)')
title(sprintf('mode %d', mode(onoffroixshift)))
subplot(2,2,4); 
histogram(onoffroiyshift, 'BinWidth', 1)
title(sprintf('mode %d', mode(onoffroiyshift)))
xlabel('Y-shift (pixels)')

% figure;
% % imshow(offline.ops.meanImg/300)
% imshow(10*offline.ops.Vcorr)
% hold on
% plot(offroictr(offiscell,2), offroictr(offiscell,1), 'r.')

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

off2ontimeinds = cat(2, off2ontimeinds, starttimeind(startfileind):endtimeind(endfileind));
end
%%이게 전 꺼

numframess2pfile = [];
for f = 1:size(suite2pfilelist,1)%여기에서 전체 tiff파일의 프레임 수를 로드하기,offline 경로
    import ScanImageTiffReader.ScanImageTiffReader.*;
    % dotindex = strfind(suite2pfilelist(f,:), '.h5');
    % fnsplit = strsplit(suite2pfilelist(f,1:dotindex-1), '/'); %이 부분은
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


off2ontimeinds = [];
for f = 1:numel(onlines2ph5folders)
starttimeind = [1; cumsum(numtimepointss2pfile(1:end-1))+1];
endtimeind = cumsum(numtimepointss2pfile);
off2ontimeinds = cat(2, off2ontimeinds, starttimeind(startid(f)):endtimeind(endid(f)));%여기에서 온라인 게 몇번째인지지 지정하는 것
end
%%여기까지
if ~( length(off2ontimeinds)==size(online.F,2) )
    error('error in matching time indices of online and offline processing')
end

onlineF = cat(1,online.F);
onoffcorrF = corr(onlineF(isneuron,:)', offline.F(:,off2ontimeinds)');
if ~isequal(size(onoffroiumdist), size(onoffcorrF))
    error('check size mismatch')
end
% sv = sort(onoffcorrF, 2, 'descend'); disp(sv(randperm(size(onoffcorrF,1),5), 1:5))
[maxon2offcorrF, imaxon2offcorrF] = max(onoffcorrF, [], 2);
[minon2offdist, iminon2offdist] = min(onoffroiumdist, [], 2);
fprintf('%.2f%% match in max F correlation and minimum distance online to offline\n', 100*mean(imaxon2offcorrF==iminon2offdist))

[maxoff2oncorrF, imaxoff2oncorrF] = max(onoffcorrF, [], 1);
[minoff2ondist, iminoff2ondist] = min(onoffroiumdist, [], 1);
fprintf('%.2f%% match in max F correlation and minimum distance offline to online\n', 100*mean(imaxoff2oncorrF==iminoff2ondist))

d = onoffroiumdist(sub2ind(size(onoffroiumdist), (1:nnz(isneuron))', imaxon2offcorrF));
figure; plot(d, maxon2offcorrF, '.')
xlabel('Distance (Pixels)')
ylabel('Max On-Off F Correlation')
title('For Each Online ROI')
xlim([0 100])

save(strcat(pathpp, 'online2offline.mat'), 'allfoldernames', 'onoffroiumdist', 'off2ontimeinds', 'onoffcorrF')
