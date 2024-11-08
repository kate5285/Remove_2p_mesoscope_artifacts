% fast way to remove stim artifcact post-suite2p
% the proper way is to remove it form the tif files itself *before* suite2p
%% Observe the data
% plot raw fluorescence trace, ordered by y-axis coordinates (leftmost strip top to bottom, second left strip top to bottom, ... , rightmost strip top to bottom)
% mark the timepoints and y-coordinates affected by the artifact with a red shade
% mark the edge zone (x-coordinates) with a deeper shade
clear all; close all; clc
% addpath('C:\Users\USER\GitHub\helperfunctions')

ses2agg = {'MU31_2/230106/','MU31_1/230106/','MU31_2/230107/','MU31_1/230107/', ...
    'MU31_1/230108/','MU31_2/230109/','MU31_1/230109/','MU31_2/230110/','MU31_1/230110/', ...
    'MU31_2/230111/','MU31_1/230111/'};
pltopt = true;

for ises = 6 % numel(ses2agg)-1:-1:7
    clearvars -except pltopt ses2agg ises
    mousedate = ses2agg{ises};
    sesclk = tic;

    mesos2ppath = ['\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts\' mousedate 'suite2p' filesep];
    load([mesos2ppath 'combined' filesep 'Fall.mat']);

    % load tif files, figure out artifact zone/pixels, apply motion correction,
    % motion correction shifts: ops.xoff, ops.yoff
    for istrip = 1:ops.nrois
        tempops = load([mesos2ppath 'plane' num2str(istrip-1) filesep 'Fall.mat'], 'ops');
        if istrip==1
            opseachstrip = tempops.ops;
        else
            opseachstrip = cat(1, opseachstrip, tempops.ops);
        end
        clearvars tempops
    end
    opseachstripxoff = cat(1,opseachstrip.xoff);
    opseachstripyoff = cat(1,opseachstrip.yoff);
    if isequal(opseachstripxoff(end,:), ops.xoff) && isequal(opseachstripyoff(end,:), ops.yoff)
    else
        warning('combined ops xoff/yoff is NOT identical to the last stripe')
    end
    if pltopt
        figure('Position', [0 100 300 300])
        subplot(1,2,1)
        histogram(opseachstripxoff)
        xlabel('X-axis motion correction')
        ylabel('# frames')
        title(mousedate, 'interpreter', 'none')
        subplot(1,2,2)
        histogram(opseachstripyoff)
        xlabel('Y-axis motion correction')
        ylabel('# frames')
        prctile( abs(opseachstripyoff(:)), 99)
        max(opseachstripyoff,[],2)
    end

    % allROImed(:,1) is vertical axis, allROImed(:,2) is horizontal axis
    allROIstat = cat(1,stat{:});
    allROImed = cat(1,allROIstat.med);

    % typical ROI size
    allROIxpix = {allROIstat.xpix};
    allROIypix = {allROIstat.ypix};
    allROINpix = cellfun(@numel, allROIxpix, 'uniformoutput', false);
    allROINpix = double(cat(1,allROINpix{:}));
    allROIxrange = cellfun(@range, allROIxpix, 'uniformoutput', false);
    allROIxrange = double(cat(1,allROIxrange{:}));
    allROIyrange = cellfun(@range, allROIypix, 'uniformoutput', false);
    allROIyrange = double(cat(1,allROIyrange{:}));
    disp([mousedate ' typical ROI size x/y/radius median/mean'])
    disp([median(allROIxrange) mean(allROIxrange);
        median(allROIyrange) mean(allROIyrange);
        median(2*sqrt(2*allROINpix/pi)) mean(2*sqrt(2*allROINpix/pi))])
    % artifact buffer zone only concered with vertical axis:
    % median(allROIyrange)/2 is 7.5 (MU31_2/230106/)
    % prctile( abs(opseachstripyoff(:)), 99) is 3 (MU31_2/230106/)
    % set buffer zone as the sum of these values +alpha
    % ybuffer = 15; % pixels
    % xbuffer = 27; % pixels



    pathpp = ['\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts_postprocessed\' mousedate];
    load([pathpp 'postprocessed.mat'], 'allfoldernames', 'exptids', 'nexpts')

    % foldername = 'SGholo';
    for ifolder = 14
        foldername = allfoldernames{14};
        % if ~( contains(foldername, 'holo') || contains(foldername, 'stimtest') )
        %     continue
        % end
        clearvars -except pltopt ses2agg ises mousedate pathpp allfoldernames exptids nexpts ifolder foldername allROImed sesclk
        disp(foldername)

        folderSIpath = ['\\shinlab\ShinLab\MesoHoloExpts\mesoholoexpts_scanimage\' mousedate foldername filesep];
        folderSItifs = dir([folderSIpath '*.tif']);
        % if numel(folderSItifs)==0
        %     continue
        % end

        load(['d:\Users\USER\Documents\MATLAB\spontaneous\', 'Fall_split.mat'])
        load(['d:\Users\USER\Documents\MATLAB\spontaneous\', 'dFFcell.mat'])
        addpath("d:\Users\USER\Documents\MATLAB\Analyze_IC_mesoholo-main\")
        import getSItifsize.*;
        fname = [folderSItifs(1).folder filesep folderSItifs(1).name];%tif가 전혀 없음
        [Lx,Ly, data] = getSItifsize(fname);

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
            if strfind(hSIh{n}{1}, 'SI.hChannels.channelSave')
                channelSave = hSIh{n}{2};
            end
        end
        if ~strcmp(channelSave, '[1;2]')
            error('violated assumption that both green and red PMT channels were recorded')
        end

        import ScanImageTiffReader.ScanImageTiffReader.*;
        reader=ScanImageTiffReader(fname);
        desc=reader.descriptions();
        fnumframes = size(desc,1); % 52 frames in first file
        temptif = reader.data();

        % following takes ~15 min for SGholo
        nonedgezone = 161:size(temptif,1)-160; % in the example session (ises=1) non-edge horizontal indices are 153-442
        vertrednonedgehorzmean = zeros(size(temptif,2),size(F,2)); % used to determine PMT artifact
        vertredhorzmax = zeros(size(temptif,2),size(F,2));
        horzredvertmax = zeros(size(temptif,1),size(F,2));
        whichtiffile = zeros(1,size(F,2));
        tcnt = 0;
        tic
        for itif = 1:numel(folderSItifs)
            fname = [folderSItifs(itif).folder filesep folderSItifs(itif).name];
            reader=ScanImageTiffReader(fname);
            desc=reader.descriptions();

            tnext = tcnt+size(desc,1)/2;
            whichtiffile(tcnt+1:tnext) = itif;

            temptif = reader.data();
            vertredhorzmax(:,tcnt+1:tnext) = squeeze(max(temptif(:,:,2:2:end),[],1));
            horzredvertmax(:,tcnt+1:tnext) = squeeze(max(temptif(:,:,2:2:end),[],2));
            vertrednonedgehorzmean(:,tcnt+1:tnext) = squeeze(mean(temptif(nonedgezone,:,2:2:end),1));

            tcnt = tnext;
            % toc
            % fprintf('%d/%d\n',itif, numel(folderSItifs))
        end
        toc
        fprintf('%d/%d\n',itif, numel(folderSItifs))
        if tcnt ~= size(F,2)
            error('number of timepoints does not match between suite2p and mesoSIpath folder')
        end


        %% Temporal interpolation & smoothing
        % 1. find edge zone and suprathreshold pixels as putative stim artifact pixels
        % 	- take motion correction into account
        % 2. where does PMT artifact happen?
        % 3. do temporal interpolation for all ROIs that include artifact pixels. maybe gaussian smoothing with std of 1 frame
        % 	- do the same thing for neuropil traces of corresponding ROIs
        %   - for simplicity, take ROIs whose centroid ('stat.med') falls inside
        %   artifact zone, add vertical buffer of 15 pixels to account for soma
        %   radius and motion correction shifts

        % what to do with PMT artifact (i.e., the darkening of lines where the stim
        % artifact appears)
        % option 1: ignore
        % option 2: NaN for all the lines where the stim artifact appears + 200-250 lines *below* artifact end
        % option 3: NaN all frames where artifact appears, align to stim *end*

        redzonethresh = 1250; % between 1250 and 2270
        redpixthresh = 400;

        ybuffer = 15; % pixels
        xbuffer = 27; % pixels

        kerwinhalf = 2; kersigma = 1;
        kergauss = normpdf( (-kerwinhalf:kerwinhalf), 0,kersigma);
        kergauss = (kergauss/sum(kergauss));

        vertredzoneframes = vertredhorzmax>redzonethresh;



        % redthresh
        % to determine whether artifact appeared in the frame or not, max over each
        % horizontal line, count nubmer of lines that are suprathrehsold.
        if pltopt
            figure('Position', [400 100 400 300])
            plot(1:size(temptif,1), horzredvertmax)
            hold on; plot(1:size(temptif,1), max(horzredvertmax,[],2), 'k-', 'Linewidth', 2)
            plot([1 size(temptif,1)], redzonethresh*[1 1], 'r--', 'Linewidth', 2)
            xlim([1 size(temptif,1)])
            xlabel('tif horz axis')
            ylabel('max red')
            title([mousedate foldername], 'interpreter', 'none')
        end

        % remap the ROI coordinates to original coordinates where all the MROIs were vertically stacked
        %{
% i.e., inverse of the following
tempdata = zeros(max(Ly), sum(Lx), size(q,3));
for istrip = 1:numel(Ly)
    tempdata(1:Ly(istrip), sum(Lx(1:istrip-1))+1:sum(Lx(1:istrip)), :) = ...
        permute(q(:,data.lines{istrip}+1,:),[2 1 3]);
end
        %}
        allROItifcoord = allROImed;
        cumLx = [0; cumsum(Lx)];
        for istrip = 1:numel(Lx)
            roisinstrip = allROImed(:,2)>=cumLx(istrip) & allROImed(:,2)<cumLx(istrip+1);
            allROItifcoord(roisinstrip,2) = allROImed(roisinstrip,2) -cumLx(istrip);
            allROItifcoord(roisinstrip,1) = allROImed(roisinstrip,1) +data.lines{istrip}(1);
        end
        iscelltifcoord = round(allROItifcoord(iscell==1,:)); % make it integer
        %{
figure;
imshow(squeeze(mean(rawSItifs(:,:,1:2:end),3)/150))
caxis([0 1])
hold on
scatter(allROItifcoord(iscell==1,1), allROItifcoord(iscell==1,2), 'rx')
colorbar
% caxis([0 1000])
% xlim([0 600])
        %}

        %% option 1: NaN only the ROIs in holo stim artifact zone (ignore PMT darkening artifact)
        edgezoneframes = true(size(horzredvertmax));
        edgezoneframes(nonedgezone,:) = false;

        iscellartzonetimepoints = false(size(dFF));
        ymax = size(vertredhorzmax,1);
        for yshift = -ybuffer:ybuffer
            ycoords = round( iscelltifcoord(:,1)+yshift ); % make it integer
            ycoords(ycoords<=1)=1;
            ycoords(ycoords>=ymax) = ymax;
            xcoords = round( iscelltifcoord(:,2) ); % make it integer

            temparttp = vertredzoneframes(ycoords,:) & edgezoneframes(xcoords,:);
            iscellartzonetimepoints = iscellartzonetimepoints | temparttp;
        end

        %{
horzredzoneframes = horzredvertmax>redzonethresh;

ycoords = iscelltifcoord(:,1) ; % make it integer
xcoords = iscelltifcoord(:,2) ; % make it integer
iscellartifacttimepoints0 = vertredzone(ycoords,:) & horzredzone(xcoords,:);

% double for loop below takes ~20 seconds
tic
iscellartifacttimepoints1 = false(size(dFF));
xmax = size(horzredvertmax,1);
ymax = size(vertredhorzmax,1);
for xshift = -xbuffer:xbuffer
    for yshift = -ybuffer:ybuffer
        ycoords = round( iscelltifcoord(:,1)+yshift ); % make it integer
        ycoords(ycoords<=1)=1;
        ycoords(ycoords>=ymax) = ymax;
        xcoords = round( iscelltifcoord(:,2)+xshift ); % make it integer
        xcoords(xcoords<=1)=1;
        xcoords(xcoords>=xmax) = xmax;

        temparttp = vertredzone(ycoords,:) & horzredzone(xcoords,:);
        iscellartifacttimepoints1 = iscellartifacttimepoints1 | temparttp;
    end
end
toc

disp([nnz(iscellartifacttimepoints0) nnz(iscellartifacttimepoints) nnz(iscellartifacttimepoints1)]) % from low to high value
all(iscellartifacttimepoints1(iscellartifacttimepoints0)) % must be true
mean(iscellartifacttimepoints(iscellartifacttimepoints0))
mean(iscellartifacttimepoints1(iscellartifacttimepoints))
        %}

        Fcell = F(iscell,:);
        Fneucell = Fneu(iscell,:); % neuropil
        spkscell = spks(iscell,:);
        Fccell = Fcell - ops.neucoeff * Fneucell;
        Fczcell = (Fccell - mean(Fccell,2))./std(Fccell,0,2);

        rawact = struct();
        rawact.Fczcell = Fczcell;
        rawact.dFF = dFF;
        rawact.spkscell = spkscell;

        % actfields = {'Fczcell', 'dFF', 'spkscell'};
        actfields = fieldnames(rawact);

        artzone = struct();
        for ii = 1:numel(actfields)

            tempnanartifact = eval(actfields{ii});
            tempnanartifact(iscellartzonetimepoints) = NaN;

            tic
            tempinterp = tempnanartifact;
            for ci = 1:nnz(iscell)
                x = tempnanartifact(ci,:);
                nanx = isnan(x);
                t    = 1:size(x,2);
                x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx), 'makima');
                tempinterp(ci,:) = x;
            end
            toc

            if size(tempinterp,2)~=size(convn(tempinterp,kergauss),2) && size(tempinterp,1)==nnz(iscell) && size(convn(tempinterp,kergauss),1)==nnz(iscell)
            else
                error('convn along wrong direction')
            end
            artzone.(actfields{ii}) = convn(tempinterp,kergauss,'same');
        end


        %% option 2: NaN for all the lines where the stim artifact appears + 200-250 lines *below* artifact
        pmtartlines = 250;
        iscellartlinetimepoints = false(size(dFF));
        ymax = size(vertredhorzmax,1);
        tic
        for yshift = -ybuffer:pmtartlines
            ycoords = round( iscelltifcoord(:,1)+yshift ); % make it integer
            ycoords(ycoords<=1)=1;
            ycoords(ycoords>=ymax) = ymax;

            temparttp = vertredzoneframes(ycoords,:);
            iscellartlinetimepoints = iscellartlinetimepoints | temparttp;
        end
        toc

        artlines = struct();
        for ii = 1:numel(actfields)

            tempnanartifact = eval(actfields{ii});
            tempnanartifact(iscellartlinetimepoints) = NaN;

            tic
            tempinterp = tempnanartifact;
            for ci = 1:nnz(iscell)
                x = tempnanartifact(ci,:);
                nanx = isnan(x);
                t    = 1:size(x,2);
                x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx), 'makima');
                tempinterp(ci,:) = x;
            end
            toc

            if size(tempinterp,2)~=size(convn(tempinterp,kergauss),2) && size(tempinterp,1)==nnz(iscell) && size(convn(tempinterp,kergauss),1)==nnz(iscell)
            else
                error('convn along wrong direction')
            end
            artlines.(actfields{ii}) = convn(tempinterp,kergauss,'same');
        end


        %% option 3: NaN all frames where artifact appears, align to stim *end*
        % frames 7-11 are artifacts. 12 onward is 'off'

        [v,c]=uniquecnt(whichtiffile);
        sv = sort(c,'descend');
        Nframespertrial = sv(2);
        numdurtrialinds = mode(c);

        redlinesperframe = sum(vertredzoneframes,1);
        if pltopt
            figure('Position', [800 100 400 300])
            hold all
            for itif = 1:numel(folderSItifs)
                plot(redlinesperframe(whichtiffile==itif))
            end
            yl = ylim;
            plot(6*[1 1], yl, 'k--')
            plot(12*[1 1], yl, 'k--')
            xlabel('frames')
            ylabel('# red lines')
            xlim([0 Nframespertrial])
            title([mousedate foldername], 'interpreter', 'none')
        end

        artstartframe = NaN(numel(folderSItifs),1);
        artendframe = NaN(numel(folderSItifs),1);
        artNlines = zeros(numel(folderSItifs),1);
        for itif = 1:numel(folderSItifs)
            artNlines(itif) = max(redlinesperframe(whichtiffile==itif));
            if nnz(redlinesperframe(whichtiffile==itif))>0 % artifact exists on this trial
                artstartframe(itif) = find(redlinesperframe(whichtiffile==itif)>0,1,'first');
                artendframe(itif) = find(redlinesperframe(whichtiffile==itif)>0,1,'last');
            end
        end

        psthfields = {'Fczcell', 'artzone_Fczcell', 'artzone_dFF', 'artzone_spkscell', 'artlines_Fczcell', 'artlines_dFF', 'artlines_spkscell'};
        psthrmart = struct();
        for f = 1:numel(psthfields)
            switch psthfields{f}
                case 'Fczcell'
                    tempF = rawact.Fczcell;
                case 'artzone_Fczcell'
                    tempF = artzone.Fczcell;
                case 'artzone_dFF'
                    tempF = artzone.dFF;
                case 'artzone_spkscell'
                    tempF = artzone.spkscell;
                case 'artlines_Fczcell'
                    tempF = artlines.Fczcell;
                case 'artlines_dFF'
                    tempF = artlines.dFF;
                case 'artlines_spkscell'
                    tempF = artlines.spkscell;
            end
            psthrmart.(psthfields{f}) = NaN(nnz(iscell), numel(folderSItifs), numdurtrialinds);
            for itif = 1:numel(folderSItifs)
                trialstartind = find(whichtiffile==itif, 1, 'first');
                temptrial = tempF(:,trialstartind:trialstartind+numdurtrialinds-1);
                psthrmart.(psthfields{f})(:,itif,:) = temptrial;
            end
        end
pathppi='d:\Users\USER\Documents\MATLAB';
        save([pathppi 'rmartifact_' foldername '.mat'], ...
            'nonedgezone', 'vertrednonedgehorzmean', 'vertredhorzmax', 'horzredvertmax', 'whichtiffile', 'allROItifcoord', ...
            'redzonethresh', 'redpixthresh', 'ybuffer', 'xbuffer', 'kergauss', 'iscelltifcoord', 'pmtartlines', ...
            'rawact', 'artzone', 'artlines', 'artstartframe', 'artendframe', 'artNlines', 'psthrmart', '-v7.3')

    end
    toc(sesclk)
end

%% reduce pmtartifact lines to ybuffer
clear all; close all; clc
ses2agg = {'MU31_2/230106/','MU31_1/230106/','MU31_2/230107/','MU31_1/230107/', ...
    'MU31_1/230108/','MU31_2/230109/','MU31_1/230109/','MU31_2/230110/','MU31_1/230110/', ...
    'MU31_2/230111/','MU31_1/230111/'};

for ises = 1:numel(ses2agg)
    clearvars -except pltopt ses2agg ises
    mousedate = ses2agg{ises};
    disp(mousedate)
    sesclk = tic;

    mesos2ppath = ['S:\MesoHoloExpts\mesoholoexpts\' mousedate 'suite2p' filesep];
    load([mesos2ppath 'combined' filesep 'Fall.mat']);

    pathpp = ['S:\MesoHoloExpts\mesoholoexpts_postprocessed\' mousedate];
    load([pathpp 'postprocessed.mat'], 'allfoldernames', 'exptids', 'nexpts')

    for ifolder = 1:numel(allfoldernames)
        foldername = allfoldernames{ifolder};
        if ~( contains(foldername, 'holo') || contains(foldername, 'stimtest') )
            continue
        end
        clearvars -except pltopt ses2agg ises mousedate pathpp allfoldernames exptids nexpts ifolder foldername allROImed sesclk
        disp(foldername)

        folderSIpath = ['S:\MesoHoloExpts\mesoholoexpts_scanimage\' mousedate foldername filesep];
        folderSItifs = dir([folderSIpath '*.tif']);
        if numel(folderSItifs)==0
            continue
        end

        load(['S:\MesoHoloExpts\mesoholoexpts\' mousedate foldername filesep 'Fall_split', nexpts{ifolder}, '.mat'])
        load(['S:\MesoHoloExpts\mesoholoexpts\' mousedate foldername filesep 'dFFcell', nexpts{ifolder}, '.mat'])

        Fcell = F(iscell,:);
        Fneucell = Fneu(iscell,:); % neuropil
        spkscell = spks(iscell,:);
        Fccell = Fcell - ops.neucoeff * Fneucell;
        Fczcell = (Fccell - mean(Fccell,2))./std(Fccell,0,2);

        load([pathpp 'rmartifact_' foldername '.mat'])
        vertredzoneframes = vertredhorzmax>redzonethresh;

        % %% option 2-2: NaN for all the lines within the stim artifact
        pmtartlines0 = ybuffer;
        iscellartlinetimepoints = false(size(dFF));
        ymax = size(vertredhorzmax,1);
        tic
        for yshift = -ybuffer:pmtartlines0
            ycoords = round( iscelltifcoord(:,1)+yshift ); % make it integer
            ycoords(ycoords<=1)=1;
            ycoords(ycoords>=ymax) = ymax;

            temparttp = vertredzoneframes(ycoords,:);
            iscellartlinetimepoints = iscellartlinetimepoints | temparttp;
        end
        toc
        actfields = fieldnames(rawact); % {'Fczcell', 'dFF', 'spkscell'}
        artlines0 = struct();
        for ii = 1:numel(actfields)

            tempnanartifact = eval(actfields{ii});
            tempnanartifact(iscellartlinetimepoints) = NaN;

            tic
            tempinterp = tempnanartifact;
            for ci = 1:nnz(iscell)
                x = tempnanartifact(ci,:);
                nanx = isnan(x);
                t    = 1:size(x,2);
                x(nanx) = interp1(t(~nanx), x(~nanx), t(nanx), 'makima');
                tempinterp(ci,:) = x;
            end
            toc

            if size(tempinterp,2)~=size(convn(tempinterp,kergauss),2) && size(tempinterp,1)==nnz(iscell) && size(convn(tempinterp,kergauss),1)==nnz(iscell)
            else
                error('convn along wrong direction')
            end
            artlines0.(actfields{ii}) = convn(tempinterp,kergauss,'same');
        end

        numdurtrialinds = size(psthrmart.Fczcell,3);
        psthfields = {'artlines0_Fczcell', 'artlines0_dFF', 'artlines0_spkscell'};
        % psthrmart = struct();
        for f = 1:numel(psthfields)
            switch psthfields{f}
                case 'artlines0_Fczcell'
                    tempF = artlines0.Fczcell;
                case 'artlines0_dFF'
                    tempF = artlines0.dFF;
                case 'artlines0_spkscell'
                    tempF = artlines0.spkscell;
            end
            psthrmart.(psthfields{f}) = NaN(nnz(iscell), numel(folderSItifs), numdurtrialinds);
            for itif = 1:numel(folderSItifs)
                trialstartind = find(whichtiffile==itif, 1, 'first');
                temptrial = tempF(:,trialstartind:trialstartind+numdurtrialinds-1);
                psthrmart.(psthfields{f})(:,itif,:) = temptrial;
            end
        end

        save([pathppi 'rmartifact_' foldername '.mat'], ...
            'nonedgezone', 'vertrednonedgehorzmean', 'vertredhorzmax', 'horzredvertmax', 'whichtiffile', 'allROItifcoord', ...
            'redzonethresh', 'redpixthresh', 'ybuffer', 'xbuffer', 'kergauss', 'iscelltifcoord', 'pmtartlines', ...
            'rawact', 'artzone', 'artlines', 'artstartframe', 'artendframe', 'artNlines', 'psthrmart', ...
            'pmtartlines0', 'artlines0', '-v7.3')
    end
    toc(sesclk)
end
