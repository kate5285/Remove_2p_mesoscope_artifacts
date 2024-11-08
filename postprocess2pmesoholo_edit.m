clear all; close all; clc

pltses = false;
sesclk = tic;

%mousedate = 'HS_Chrome2f_3/220923/';
%mousedate = 'MU31_2/230106/'; % drew V1 border
%mousedate = 'MU31_1/230106/'; % drew V1 border
%mousedate = 'MU31_2/230107/'; % drew V1 border
%mousedate = 'MU31_1/230107/'; % drew V1 border
%mousedate = 'MU31_1/230108/'; % drew V1 border
mousedate = 'MU31_2/230109/'; % drew V1 border
%mousedate = 'MU31_1/230109/'; % drew V1 border
%mousedate = 'MU31_2/230110/'; % drew V1 border
%mousedate = 'MU31_1/230110/'; % drew V1 border
%mousedate = 'MU31_2/230111/'; % drew V1 border
% mousedate = 'MU31_1/230111/'; % drew V1 border

disp(mousedate)
drivepath = '\\shinlab\ShinLab\MesoHoloExpts\';
mesoSIpath = [drivepath 'mesoholoexpts_scanimage/' mousedate];
onlinepath = [drivepath 'mesoholoexpts_scanimage/' mousedate 'ClosedLoop_justgreen/'];
pathvisstim = [drivepath 'mesoholoexpts_scanimage/' mousedate 'visstiminfo/'];
path2p = [drivepath 'mesoholoexpts/' mousedate];
pathpp = [drivepath 'mesoholoexpts_postprocessed/' mousedate];

%load(strcat(pathpp, 'offline_params.mat'))

if exist(strcat(path2p, 'suite2p/combined'), 'dir') % multi plane imaging
    offline = load(strcat(path2p, 'suite2p/combined/Fall.mat')); % F, Fneu, iscell, ops, redcell, spks, stat
else % single plane imaging
    offline = load(strcat(path2p, 'suite2p/plane0/Fall.mat'));
end
load([pathpp 'offline_params.mat'])

offstat = cell2mat(offline.stat);
offmed = double( cat(1, offstat.med) );

fprintf('%d/%d ROIs iscell\n', nnz(offline.iscell(:,1)==1), size(offline.iscell,1))

%% save F, spks and dF/F for each exptid. Brace yourself, this takes a while!
% dF/F. To calculate the dF/F for each fluorescence trace,  we first calculated 
% baseline fluorescence by using a median filter of width 5,401 samples (180 s). 
% We then calculated the change in fluorescence relative to baseline fluorescence (?F), 
% divided by baseline fluorescence (F). 
% To prevent very small or negative baseline fluorescence, we set the 
% baseline as the maximum of the median filter-estimated baseline and the 
% s.d. of the estimated noise of the fluorescence trace.

% frame rate for mesoscope is too slow for mode -- use median instead

offstat = cell2mat(offline.stat);
offmed = double( cat(1, offstat.med) );

ops = offline.ops;
stat = offline.stat;

iscell = boolean(offline.iscell(:,1));
redcell = boolean(offline.redcell(:,1));

Nneurons = nnz(iscell);
notredneurons = find(redcell(iscell)==0);
redneurons = find(redcell(iscell)==1);

% after manual curation
pathppi='d:\Users\USER\Documents\MATLAB\';
save([pathppi 'postsuite2p_params.mat'], 'offmed', 'iscell', 'redcell', 'redneurons')

for ii=14
    if ii==0
        tempexptid = 'spontaneous';
        tempnexp = '';
    else
        tempexptid = exptids{ii};
        tempnexp = nexpts{ii};
    end
    if ~strcmp(tempexptid, 'spontaneous')
        exptidn = strcat(exptids{ii}, '_', nexpts{ii});
        if contains(exptids{ii}, 'holo') || contains(allfoldernames{ii}, 'stimtest') || contains(exptids{ii}, 'ref')
        elseif ~isfield(vis, exptidn)
            warning('%s does not have a vis file', exptidn)
        elseif vis.(exptidn).numtrials+1 ~= nnz(whichfolder==ii)
            warning('%s number of Scanimage tif files (%d) should match number of trials plus one (%d+1)', ...
                exptidn, nnz(whichfolder==ii), vis.(exptidn).numtrials)
        end
    end
    
    if nnz(whichfolder==ii)==0
        continue
    end
    
ntpf = find(whichfolder==ii);
sesstartind = sum(numtimepointss2pfile(1:ntpf(1)-1))+1;
sesendind = sum(numtimepointss2pfile(1:ntpf(end)));
sesframeinds = sesstartind:sesendind;
    
    % split into each exptid
    F = offline.F(:, sesframeinds);
    Fneu = offline.Fneu(:, sesframeinds);
    spks = offline.spks(:, sesframeinds);
    if ~exist(strcat(path2p, tempexptid), 'dir')
        mkdir(strcat(path2p, tempexptid))
    end
    save(strcat(pathppi, tempexptid, '/Fall_split', tempnexp, '.mat'), 'F', 'Fneu', 'iscell', 'ops', 'redcell', 'spks', 'stat')
    
    % dF/F, z-scored dF/F
    Fcell = F(iscell,:);
    Fneucell = Fneu(iscell,:); % neuropil
    spkscell = spks(iscell,:);
    
    Fccell = Fcell - ops.neucoeff * Fneucell;
    % Fzcell = (Fccell-mean(Fccell,2))./std(Fccell,0,2);
    numbaseframes = min(round(ops.fs*180), size(Fccell,2)); % 180seconds; or the entire session if session length is less than 3 min (this happens in subretinotopy)
    
    
    % the following method will be used starting September 2020
    % only do rolling average if there's a significant difference in the
    % first 3 min vs last 3 min
%     p = signrank(mean(Fcell(:,1:numbaseframes),2), mean(Fcell(:,end-numbaseframes+1:end),2));
    % turns  out signrank is too sensitive, so using a different measure
    if round(ops.fs*180)>size(Fccell,2)
        % 180seconds; or the entire session if session length is less than 3 min (this happens in subretinotopy)
        rolling = false;
    else
        Fbasefirst = median(Fccell(:,1:numbaseframes), 2);
        Fbaselast = median(Fccell(:,end-numbaseframes+1:end), 2);
        prctgt = 100*mean(Fbasefirst>Fbaselast);
        fprintf('Fbasefirst vs Fbaselast signrank p=%.4f\n    %.0f%% got dimmer\n', signrank(Fbasefirst, Fbaselast), prctgt)
        rolling = prctgt>95||prctgt<5;
    end
    
    if rolling
        % do rolling mode if greater than 95% of the cells or less than
        % 5% of the cells have higher fluorescence in the first 3 min
        % compared to last 3 min
        disp('rolling mode: this will take a while')
        
        % rolling window binned at every timepoint
        medindrow1= [zeros(1,ceil(numbaseframes/2)) 1:size(Fccell,2)-numbaseframes (size(Fccell,2)-numbaseframes)*ones(1,floor(numbaseframes/2))];
        if length(medindrow1) ~= size(Fccell, 2)
            error('check medindrow1')
        end
        medind = [1:numbaseframes]' + medindrow1;
        Fbase=zeros(size(Fccell));
        tic
        for ci=1:Nneurons
            tempF = Fccell(ci,:);
            tempFbase = median(tempF(medind), 1);
            Fbase(ci,:) = tempFbase;
        end
        toc % takes ~3 times longer than doing median        
        
        % rolling window binned at 1/2 binwidth
%         medbinoverlap = 2;
%         numbaseframes = round(numbaseframes/medbinoverlap)*medbinoverlap;  % make it an even number
%         medind = [1:numbaseframes]' + (0:numbaseframes/medbinoverlap:size(Fccell,2)-numbaseframes);
%         Fbase=zeros(size(Fccell));
%         tic
%         for ci=1:Nneurons
%             tempF = Fccell(ci,:);
%             tempFbase = median(tempF(medind), 1);
%             Fbase(ci,1:length(tempFbase)*numbaseframes/medbinoverlap) = ...
%                 reshape(repmat(tempFbase,numbaseframes/medbinoverlap,1),[],1)';
%             Fbase(ci,length(tempFbase)*numbaseframes/medbinoverlap+1:end) = tempFbase(end);
%         end
%         toc % much much faster                
    else
%         Fbase = mode(bw*(-0.5+round(Fccell/bw +0.5)), 2);
        Fbase = median(Fccell, 2);
        Fbase = repmat(Fbase, 1,size(Fccell,2));
    end    
    if ~isequal(size(Fccell), size(Fbase))
        error('F0 and F must be same size. check code')
    end
    dFF = (Fccell-Fbase)./Fbase;
    
    save(strcat(pathppi, tempexptid, '/dFFcell', tempnexp, '.mat'), 'dFF', 'spkscell', 'notredneurons', 'redneurons');
end
clear('offline', 'F', 'Fneu','spks')



%% psth for each exptid. this takes less than 10s
% load(strcat(path2p, 'suite2p/combined/Fall.mat'), 'ops', 'iscell')
% Nneurons = nnz(iscell);

psthall = struct();
Rall = struct(); % response during stimulus duration
validexpts = false(numel(allfoldernames), 1);
tic
for ii= 14
    try
        load(strcat(path2p, exptids{ii}, '/dFFcell', nexpts{ii}, '.mat')) % 'dFF', 'dFFz', 'spkscell', 'notredneurons', 'redneurons'
    catch
        continue
    end
    
    load(strcat(path2p, exptids{ii}, '/Fall_split', nexpts{ii}, '.mat'))
    Fcell = F(iscell,:);
    % Fczcell = (Fccell - mean(Fccell,2))./std(Fccell,0,2);

    Fneucell = Fneu(iscell,:); % neuropil
    if ops.neucoeff ~= 0.7
        warning('ops.neucoeff is not 0.7?!')
    end
    Fccell = Fcell - ops.neucoeff * Fneucell;
    Fczcell = (Fccell - mean(Fccell,2))./std(Fccell,0,2);

    exptidn = strcat(exptids{ii}, '_', nexpts{ii});
    disp(exptidn)
    
    Rall.(exptidn).Fcz_sesavg = mean(double(Fczcell),2);
    Rall.(exptidn).dFF_sesavg = mean(double(dFF),2);
    Rall.(exptidn).spks_sesavg = mean(double(spkscell),2);    
    
    Rall.(exptidn).Fcz_sesstd = std(double(Fczcell),0,2);
    Rall.(exptidn).dFF_sesstd = std(double(dFF),0,2);
    Rall.(exptidn).spks_sesstd = std(double(spkscell),0,2);
        
    if strcmp(exptids{ii}, 'blankscreen') || contains(exptids{ii}, 'spon') || contains(exptids{ii}, 'ref')
        %|| isempty(str2num(nexpts{ii})) %|| strcmp(exptids{ii}, 'RFcircleCI0')
        psthall.(exptidn).Fcz = Fczcell;
        psthall.(exptidn).dFF = dFF;
        psthall.(exptidn).spks = spkscell;
        continue
    end    
    
    % note, the last frame of each file is trialtriginds
    trialtriginds = cumsum(numtimepointss2pfile(whichfolder==ii));
    % the first few files might be mis-triggers from SI...
    while trialtriginds(1)==1
        trialtriginds(1)=[];
    end
%     if ~contains(foldername{ii}, 'holo')
        if trialtriginds(end) == size(dFF,2)
            trialtriginds(end) = [];
        else
            error('number of frames mismatch')
        end
%     end
%     if contains(exptids{ii}, 'retinotopy')
%         trialtrigframes = trialtrigframes(2:end-10);
%     end
    
    if numel(trialtriginds)==0
        warning(strcat(exptidn, ' had no triggers registered'))
        continue
    end
    
    %%%%%%%% RESUME HERE: INHERIT HOLO TRIAL STRUCTURE %%%%%%%%%%%%%%%
    if ~isfield(vis, exptidn) || contains(exptids{ii}, 'holo') || contains(allfoldernames{ii}, 'stimtest')
        numpreinds = 0;
        numpostinds = 0;
        numdurtrialinds = mode(numtimepointss2pfile(whichfolder==ii));
        fprintf('numdurtrialinds %d\n', numdurtrialinds)
        durstiminds = 1:numdurtrialinds; %%%%%%%%%%%%%%%%%%%%%
        %warning('THIS PART NEEDS EDITING')
        
        % in my vis stim, there is a spontaneous period in the beginning --
        % hence the first trial is the second file
        % in holography, the first trial is the first file
        trialstartind = [1; trialtriginds];
        numrectrials = numel(trialstartind); % number of recorded trials
    else
        if ~isequal(numel(trialtriginds), vis.(exptidn).numtrials)
            warning('%s number of 2p triggers (%d) should match number of trials (%d)', ...
                exptidn, numel(trialtriginds), vis.(exptidn).numtrials)
        end
        
        % PSTH
        numpreinds = floor(ops.fs * 2 );
        numpostinds = floor(ops.fs * 2 );
        numdurtrialinds = floor(ops.fs * vis.(exptidn).durvisstim);
        
        durstiminds = numpreinds+2:numpreinds+numdurtrialinds+1;
        trialstartind = trialtriginds;
        numrectrials = numel(trialstartind); % number of recorded trials
    end
    
    %trialtimeline = -(numpreinds + numdurtrialinds):numpostinds;
    trialtimeline = -numpreinds:numdurtrialinds+numpostinds;
    
    %if numrectrials<vis.(exptidn).numtrials
    while trialstartind(end)+trialtimeline(end) > size(dFF,2)
        numrectrials = numrectrials-1;
        trialstartind = trialstartind(1:numrectrials);
    end
    %end
        
    psthall.(exptidn).trialstartind = trialstartind;
        
    psthtrialinds = trialstartind + trialtimeline;
    psthtimeline = 1/ops.fs * [-numpreinds:numdurtrialinds+numpostinds];
    
    % psthall.(exptid).numpreinds = numpreinds;
    % psthall.(exptid).numpostinds = numpostinds;
    psthall.(exptidn).numdurtrialinds = numdurtrialinds;
    psthall.(exptidn).psthtimeline = psthtimeline;
    
    psthall.(exptidn).Fcz = zeros([Nneurons size(psthtrialinds)]);
    psthall.(exptidn).dFF = zeros([Nneurons size(psthtrialinds)]);
%     psthall.(exptidn).dFFz = zeros([Nneurons size(psthtrialinds)]);
    psthall.(exptidn).spks = zeros([Nneurons size(psthtrialinds)]);
    for ci = 1:Nneurons % <3 seconds for 1000cells
        tempFcz = Fczcell(ci,:);
        psthall.(exptidn).Fcz(ci,:,:) = tempFcz(psthtrialinds);
        
        tempdFF = dFF(ci,:);
        psthall.(exptidn).dFF(ci,:,:) = tempdFF(psthtrialinds);
        
%         tempdFFz = dFFz(ci,:);
%         psthall.(exptidn).dFFz(ci,:,:) = tempdFFz(psthtrialinds);

        tempspks = spkscell(ci,:);
        psthall.(exptidn).spks(ci,:,:) = tempspks(psthtrialinds);
    end
    

    if isfield(vis, exptidn)
        Rall.(exptidn).Fcz_sesbegin = mean(Fczcell(:,1:trialtriginds(1)-1),2);
        Rall.(exptidn).dFF_sesbegin = mean(dFF(:,1:trialtriginds(1)-1),2);
        Rall.(exptidn).spks_sesbegin = mean(spkscell(:,1:trialtriginds(1)-1),2);
        
        Rall.(exptidn).Fcz_sesend = mean(Fczcell(:,trialtriginds(end)+1:end),2);
        Rall.(exptidn).dFF_sesend = mean(dFF(:,trialtriginds(end)+1:end),2);
        Rall.(exptidn).spks_sesend = mean(spkscell(:,trialtriginds(end)+1:end),2);
    end
    
    %     Rall.(exptidn).dFF = squeeze(mean(psthall.(exptidn).dFF(:,:,numpreinds+1:numpreinds+floor(ops.fs * 2)+1),3)); % 2 seconds after stim onset
    % %     Rall.(exptidn).dFFz = squeeze(mean(psthall.(exptidn).dFFz(:,:,numpreinds+1:numpreinds+floor(ops.fs * 2)+1),3)); % 2 seconds after stim onset
    %     Rall.(exptidn).spks = squeeze(mean(psthall.(exptidn).spks(:,:,numpreinds+1:numpreinds+floor(ops.fs * 2)+1),3)); % 2 seconds after stim onset
    for wf = 0:2
        switch wf
            case 0
                whichF = 'Fcz';
            case 1
                whichF = 'dFF';
            case 2
                whichF = 'spks';
        end
        % stim duration response
%         Rall.(exptidn).(whichF) = squeeze(mean(psthall.(exptidn).(whichF)(:,:,numpreinds+1:numpreinds+numdurtrialinds ),3));        
        Rall.(exptidn).(whichF) = squeeze(mean(psthall.(exptidn).(whichF)(:,:,durstiminds ),3));        
    end
    
    if contains(allfoldernames{ii}, 'holo')
        %warning('IMPLEMENT VALIDEXPTS FOR HOLO?')
    elseif ~isfield(vis, exptidn)
            validexpts(ii) = false;
    else
        if numrectrials >= vis.(exptidn).numtrials*0.8
            validexpts(ii) = true;
        else
            validexpts(ii) = false;
        end
    end
end
toc

save(strcat(pathppi, 'postprocessed.mat'), 'allfoldernames', ...
    'validexpts', 'exptids', 'nexpts', 'vis', 'Rall', '-v7.3')
save(strcat(pathppi, 'postprocessed_psth.mat'), 'psthall', '-v7.3')

toc(sesclk)

% [Roff, Rbs, R1s, Rmax, Rrange] = computeSensoryResponses(psthall);

%%
if pltses 
figure
for ii= 14    
    exptidn = strcat(exptids{ii}, '_', nexpts{ii});
    disp(exptidn)
    
    if ~exist(strcat(path2p, exptids{ii}, '/dFFcell', nexpts{ii}, '.mat'), 'file')
        continue
    end
    
    if strcmp(exptids{ii}, 'blankscreen') %|| isempty(str2num(nexpts{ii})) %|| strcmp(exptids{ii}, 'RFcircleCI0')
        continue
    end
    
    if ~isfield(vis, exptidn)
        warning(strcat('skipped ', exptidn))
        continue
    end
    
    
    trialtriginds = cumsum(numtimepointss2pfile(whichfolder==ii));    
    trialtriginds(end) = [];
    
    numrectrials = numel(trialtriginds); % number of recorded trials
    
    if numel(trialtriginds)==0
        warning(strcat(exptidn, ' had no triggers registered'))
        continue
    end
    
    if nnz(diff(trialtriginds)<0) % frame number is mod 2^16
        resetinds = [find(diff(trialtriginds)<0); numel(trialtriginds)];
        for ri = 1:nnz(diff(trialtriginds)<0)
            trialtriginds(resetinds(ri)+1:resetinds(ri+1)) = ri*2^16 + trialtriginds(resetinds(ri)+1:resetinds(ri+1));
        end
    end
    if ~isequal(numel(trialtriginds), vis.(exptidn).numtrials)
        warning('number of 2p triggers must be same as the number of trials')
    end
    
    subplot(4,4,ii)
    if isfield(vis.(exptidn), 'Tstartvisstim')
        numend = min([numrectrials vis.(exptidn).numtrials]);
    plot(vis.(exptidn).Tstartvisstim(1:numend), trialtriginds(1:numend), '.')
    end
    title(exptidn, 'interpreter', 'none')
end
end

%%
%{
% load(strcat(pathpp, 'postprocessed.mat'))
% load(strcat(pathpp, 'postprocessed_psth.mat'))
% [Roff, Rbs, R1s, Rmax, Rrange] = computeSensoryResponses(psthall);

disp('analyzeRFcircleCI')
tic
analyzeRFcircleCI(pathpp, exptids, nexpts, validexpts, vis, Rall, ...
    Roff, Rbs, R1s, Rmax, Rrange, psthall, false)
toc

disp('analyzeGratings')
tic
analyzeGratings(pathpp, exptids, nexpts, validexpts, vis, Rall, ...
    Roff, Rbs, R1s, Rmax, Rrange, psthall);
toc

disp('analyzeStaticICtx')
tic
analyzeStaticICtx(pathpp, exptids, nexpts, validexpts, vis, Rall, ...
    Roff, Rbs, R1s, Rmax, Rrange, psthall); %, psthallus);
toc
%}
%%
% addpath(strcat(drivepath, 'CODE/Analyze_IC_mesoscope'))
% 
% analyze_retinotopy_meso
% plotFOV_retinotopy
% check_retinotopy_pixel_vs_ROI 
% analyze_visfieldsign
% 
% plotFOV_SPgest
% plotFOV_RFcircleCI

