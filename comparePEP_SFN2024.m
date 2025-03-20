%% Script to generate outputs for Figure 1, adapted from code for SFN 2024 abstract (E. Bartoli)
%
%    If this code is used in a publication, please cite the manuscript:
%    "Proper reference selection and re-referencing to mitigate bias in single pulse electrical stimulation data"
%    by H. Huang, J.A. Adkinson, M.A. Jensen, M. Hasen, I.A. Danstrom, K.R. Bijanki, N.M. Gregg, K.J. Miller,
%    S.A. Sheth, D. Hermes, E. Bartoli.
%
%    comparePEP_SFN2024.m
%    Copyright (C) 2025 Harvey Huang
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Configured to run separately for Baylor and Mayo data, depending on global "dataStr" variable
% Loads SPES data, applies CARLA, calculates differences between the 2 references before and after CARLA
% Calculates CRP outputs for subject 2

%% 1) Load data.
% dataStr configured in main script (sub1 = Baylor, sub2 = Mayo)

if strcmpi(dataStr, 'baylor')

    load('0023_PEP_RT2aHa_1.mat')
    
    % hipp reference
    DataW = Data; clear Data
    
    % dataset 2: subgal reference, stim to same hippocampal electrode
    load('0074_PEP-SubgalealRef_RT2aHa_1.mat')
    
    % subgaleal reference
    DataG = Data; clear Data
    
    % load also the colormap that we will use for correlations:
    
    %load('spect_cmap.mat')
    
    % remove bad chans
    
    srate = 2000;
    tt = xAxis/1000; % units of seconds
    
    all = 1:size(DataW,1);
    bad = [6, 7, 85, 109, 110, 122]; % bad channels
    good_mask = ~ismember(all,bad)';
    
    % keep track of which are the good channels
    ch_ind = all(good_mask);
    ch_names = MontageInfo.Current.Label(good_mask);
    nChsGood = length(ch_ind);
    
    % pull out good channels
    DataWGood = DataW(good_mask, :, :);
    DataGGood = DataG(good_mask, :, :);

elseif strcmpi(dataStr, 'mayo_LK1-LK2')
    
    stim_site = 'LK1-LK2';
    badTrs1 = [8, 9]; % bad trials for Data1 and Data2
    badTrs2 = [];
    
    load(sprintf('sub-mayo_dataReferencing_%s.mat', stim_site)); % srate is loaded, 4800
        
    % rename vars to match baylor data.
    bad = find(~strcmp(channelsInfo.status, 'good') | ~strcmp(channelsInfo.type, 'SEEG'));
    
    % exclude stimulated & stim-adjacent channels
    stim_chs = getNeighborChs(split(stim_site, '-'), 1);
    bad = [bad; find(ismember(channelsInfo.name, stim_chs))];
    
    % other manually annotated bad channels. 1:45 are all LA-LC, 65 and 97 are the references for 2nd and then 1st datasets
    bad = unique([bad; (1:45)'; 65; 97; 196]);
    
    good_mask = ~ismember(1:size(data1, 1), bad)';
    good_mask_num = find(good_mask);
    
    DataWGood = data2(good_mask, :, :); % contaminated
    DataGGood = data1(good_mask, :, :); % more neutral ref
    nChsGood = size(DataWGood, 1);
    
    % remove badtrials
    DataWGood(:, :, badTrs2) = [];
    DataGGood(:, :, badTrs1) = [];
    
    tt = tt1;
    xAxis = tt*1000; % ms scale, to match baylor data

end

%% 2) run CARLA on both datasets --------------------------------------------------------------
% Note this is non-deterministic and final order, chsUsed may differ on iterations

% CARLA on contaminated reference
[VoutW, CARW, statsW] = ccep_CARLA_modified(tt, DataWGood, srate);

% CARLA on good reference
[VoutG, CARG, statsG] = ccep_CARLA_modified(tt, DataGGood, srate);

% save for faster loading
save(fullfile('output', sprintf('carla_output_%s.mat', dataStr)), 'VoutW', 'VoutG', 'CARW', 'CARG', 'statsW', 'statsG');

%% 3) save intermediate trial-avged output to use for next step Latency_FindPeaks_EB.m

mat1 = mean(DataWGood, 3)';
mat2 = mean(DataGGood, 3)';
mat3 = mean(VoutW, 3)';
mat4 = mean(VoutG, 3)';
time = tt*1e3; % units of ms
save(fullfile('output', sprintf('intermediate_%s.mat', dataStr)), 'mat1', 'mat2', 'mat3', 'mat4', 'time', 'dataStr');

% PEP visualization - Figure 1 Top

xlims = [-100, 300];

figure('Position', [200, 200, 1000, 500]);
subplot(141)
imagesc(time,[1:size(mat1,2)],mat1');clim([-150 150]);xlim(xlims);title('WM Ref'); hold on
plot([0 0],[1 size(mat1,2)],'Color',[0 0 0])
subplot(142)
imagesc(time,[1:size(mat1,2)],mat2');clim([-150 150]);xlim(xlims);title('Neutral Ref'); hold on
plot([0 0],[1 size(mat1,2)],'Color',[0 0 0]); yticks([]);
subplot(143)
imagesc(time,[1:size(mat1,2)],mat3');clim([-150 150]);xlim(xlims);title('WM-Rereferenced'); hold on
plot([0 0],[1 size(mat1,2)],'Color',[0 0 0]); yticks([]);
subplot(144)
imagesc(time,[1:size(mat1,2)],mat4');clim([-150 150]);xlim(xlims);title('Neut-Rereferenced'); hold on
plot([0 0],[1 size(mat1,2)],'Color',[0 0 0]); yticks([]);

saveas(gcf, fullfile('output', sprintf('meanV_CARLA_horz_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('meanV_CARLA_horz_%s', dataStr)), 'svg');

% Save a colorbar

figure('Position', [400, 400, 300, 200]);
imagesc(time,[1:size(mat1,2)],mat4');clim([-150 150]);xlim(xlims);title('Neut-Rereferenced'); hold on
plot([0 0],[1 size(mat1,2)],'Color',[0 0 0]); yticks([]);
colorbar;
saveas(gcf, fullfile('output', 'meanV_colorbar'), 'svg');

%% 4) Plot the CAR for each reference & calculate RMS, for Figure 2

if strcmp(dataStr, 'baylor')
    ylims = [-600, 600];
    boxlims = [0, 250]; boxticks = 0:50:250;
elseif strcmp(dataStr, 'mayo_LK1-LK2')
    ylims = [-400, 400];
    boxlims = [0, 100]; boxticks = 0:20:100;
end

figure('Position', [200, 200, 400, 600]);
subplot(2, 1, 1); hold on; % for the contaminated reference
plot(tt, CARW, 'Color', [0.5, 0.5, 0.5]);
plot(tt, mean(CARW, 2), 'k-', 'LineWidth', 1);
xlim([-0.1, 0.3]); ylim(ylims); xticks(-0.1:0.1:0.3);
subplot(2, 1, 2); hold on; % for the neutral reference
plot(tt, CARG, 'Color', [0.5, 0.5, 0.5]);
plot(tt, mean(CARG, 2), 'k-', 'LineWidth', 1);
xlim([-0.1, 0.3]); ylim(ylims); xticks(-0.1:0.1:0.3);
saveas(gcf, fullfile('output', sprintf('CA_compare_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('CA_compare_%s', dataStr)), 'svg');

% Calculate RMS of the mean signal
rmsWmean = rms(mean(CARW(xAxis >= 9 & xAxis < 200, :), 2));
rmsGmean = rms(mean(CARG(xAxis >= 9 & xAxis < 200, :), 2));
fprintf('RMS of mean: %0.1f, %0.1f\n', rmsWmean, rmsGmean);

%% 5) T-test heatmap to find time points different between datasets

cmkjm = ieeg_kjmloccolormap();

tstatsRaw = nan(nChsGood, length(tt));
tstatsCarla = nan(nChsGood, length(tt));
for ii = 1:nChsGood
    fprintf('.');
    [~, ~, ~, statsRaw] = ttest2(squeeze(DataWGood(ii, :, :))', squeeze(DataGGood(ii, :, :))');
    tstatsRaw(ii, :) = statsRaw.tstat;
    [~, ~, ~, statsCarla] = ttest2(squeeze(VoutW(ii, :, :))', squeeze(VoutG(ii, :, :))');
    tstatsCarla(ii, :) = statsCarla.tstat;
end
fprintf('\n');

% plot heatmaps of tstats
clims = [-10, 10];
xlims = [-0.1, 0.3]; % time points to match baylor data
figure('Position', [200, 200, 500, 500]);
subplot(1, 2, 1);
imagesc(tt, 1:nChsGood, tstatsRaw); xlim(xlims); clim(clims); title('Raw-Raw'); colormap(cmkjm);
xline(0, 'Color', 'k');
subplot(1, 2, 2);
imagesc(tt, 1:nChsGood, tstatsCarla); xlim(xlims); clim(clims); title('Carla-Carla'); colormap(cmkjm);
xline(0, 'Color', 'k');
yticks([]);
saveas(gcf, fullfile('output', sprintf('tstats_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('tstats_%s', dataStr)), 'svg');

% save colorbar with separate plot
figure('Position', [400, 400, 300, 200]);
imagesc(tt, 1:nChsGood, tstatsCarla); xlim(xlims); clim(clims); title('for colorbar'); colormap(cmkjm);
colorbar;
saveas(gcf, fullfile('output', 'tstats_colorbar'), 'svg');


%% 6) Boxplots of mean abs(tstat)

meanLims = [0.009, 0.2];

% Quantify sum of abs(tstats) using boxplots
tstatsMeanRaw = mean(abs(tstatsRaw(:, tt >= meanLims(1) & tt < meanLims(2))), 2);
tstatsMeanCarla = mean(abs(tstatsCarla(:, meanLims(1) & tt < meanLims(2))), 2);
tstatsMeanX1 = mean(abs(tstatsX1(:, meanLims(1) & tt < meanLims(2))), 2);
tstatsMeanX2 = mean(abs(tstatsX2(:, meanLims(1) & tt < meanLims(2))), 2);

figure('Position', [200, 200, 150, 300]);
boxplot([tstatsMeanRaw, tstatsMeanCarla], 'Symbol', 'r.', 'Widths', 0.5);
if strcmpi(dataStr, 'baylor') % dumb y-lim switching
    ylim([0, 5]); yticks(1:5);
elseif strcmpi(dataStr, 'mayo_LK1-LK2')
    ylim([0, 2]); yticks(0.5:0.5:2);
end
xticklabels({'Raw', 'Carla'});
saveas(gcf, fullfile('output', sprintf('tstats_mean_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('tstats_mean_%s', dataStr)), 'svg');

% signed rank test that paired difference in abs(tstats) decreased with referencing
p = signrank(tstatsMeanCarla, tstatsMeanRaw, 'tail', 'left');

if strcmpi(dataStr, 'baylor'), return; end

%% 7) CRP to find significant connections (Sub 2 only)

% CRP duration for each data set
dursW = nan(nChsGood, 1);
dursG = nan(nChsGood, 1);
dursWCarla = nan(nChsGood, 1);
dursGCarla = nan(nChsGood, 1);

% CRP start time and end time
t1 = 0.009;
t2 = 1;

for ii = 1:nChsGood
    fprintf('.');
    % CRP on contam dataset, without re-ref
    [crp_parms, crp_projs] = CRP_method(squeeze(DataWGood(ii, tt >= t1 & tt < t2, :)), tt(tt >= t1 & tt < t2));
    max_tpts = length(crp_projs.proj_tpts); % number of downsampled time points
    
    h = ttest(crp_parms.al); % if coefficients are significant and duration isnt max, tR is meaningful
    if h && crp_projs.tR_index ~= max_tpts, dursW(ii) = crp_parms.tR; end

    % do same with other datasets
    [crp_parms, crp_projs] = CRP_method(squeeze(DataGGood(ii, tt >= t1 & tt < t2, :)), tt(tt >= t1 & tt < t2));
    h = ttest(crp_parms.al); % if coefficients are significant, store duration
    if h && crp_projs.tR_index ~= max_tpts, dursG(ii) = crp_parms.tR; end
    [crp_parms, crp_projs] = CRP_method(squeeze(VoutW(ii, tt >= t1 & tt < t2, :)), tt(tt >= t1 & tt < t2));
    h = ttest(crp_parms.al); % if coefficients are significant, store duration
    if h && crp_projs.tR_index ~= max_tpts, dursWCarla(ii) = crp_parms.tR; end
    [crp_parms, crp_projs] = CRP_method(squeeze(VoutG(ii, tt >= t1 & tt < t2, :)), tt(tt >= t1 & tt < t2));
    h = ttest(crp_parms.al); % if coefficients are significant, store duration
    if h && crp_projs.tR_index ~= max_tpts, dursGCarla(ii) = crp_parms.tR; end
end
fprintf('\n');

% this normalization currently considers nan as entries
figure('Position', [200, 200, 600, 400]); hold on
histogram(dursW, 0:0.05:t2, 'Normalization', "probability", 'FaceColor', [0.78, 0.71, 0.3]);
histogram(dursG, 0:0.05:t2, 'Normalization', "probability", 'FaceColor', [0.37, 0.76, 0.83]);
histogram(dursWCarla, 0:0.05:t2, 'Normalization', "probability", 'FaceColor', [0.49, 0.09, 0]);
histogram(dursGCarla, 0:0.05:t2, 'Normalization', "probability", 'FaceColor', [0, 0.19, 0.60]);
ylim([0, 0.8]);

saveas(gcf, fullfile('output', sprintf('CRP_durations_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('CRP_durations_%s', dataStr)), 'svg');

save(fullfile('output', sprintf('CRP_durations_%s.mat', dataStr)), 'dursW', 'dursG', 'dursWCarla', 'dursGCarla');
