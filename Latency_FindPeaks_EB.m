%% Script to generate outputs for Figure 2
%
%    If this code is used in a publication, please cite the manuscript:
%    "Proper reference selection and re-referencing to mitigate bias in single pulse electrical stimulation data"
%    by H. Huang, J.A. Adkinson, M.A. Jensen, M. Hasen, I.A. Danstrom, K.R. Bijanki, N.M. Gregg, K.J. Miller,
%    S.A. Sheth, D. Hermes, E. Bartoli.
%
%    Latency_FindPeaks_EB.m
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
% Reads intermediates from comparePEP_SFN2024, compute peaks, calculate correlations between channels
% Calculate metrics like means, entropy of distribution, mode size

%% 1) Load Mayo, Baylor Datasets intermediates
% dataStr configured in main script

% mat1-4 (3, 4 are re-ref)
ndatasets = 4;
if strcmpi(dataStr, 'baylor')
    load('intermediate_baylor.mat');
elseif strcmpi(dataStr, 'mayo_LK1-LK2')
    load('intermediate_mayo_LK1-LK2.mat');
end

% load the spect cmap
load('spect_cmap.mat')

%% 2) Correlation between channels
% 9-200 ms to avoid artifact and match window used for t-stat heatmaps

timepost=time>=9.000 & time<=200;

[r1,p1]=corr(mat1(timepost,:),'type','Spearman');
[r2,p2]=corr(mat2(timepost,:),'type','Spearman');
[r3,p3]=corr(mat3(timepost,:),'type','Spearman');
[r4,p4]=corr(mat4(timepost,:),'type','Spearman');

% extract upper triangular values (excluding diagonal) and get mean correlation
m  = triu(true(size(r1)),1);
fprintf('Mean rhos = %0.2f, %0.2f, %0.2f, %0.2f\n', mean(r1(m)), mean(r2(m)), mean(r3(m)), mean(r4(m)));

% how many pairs are sig (percentage of total pairs) for original data
tot_pairs = sum(sum(tril(r1,-1)~=0));
frac1 = 100*sum(sum(tril(p1<0.01,-1)))/tot_pairs;
frac2 = 100*sum(sum(tril(p2<0.01,-1)))/tot_pairs;
frac3 = 100*sum(sum(tril(p3<0.01,-1)))/tot_pairs;
frac4 = 100*sum(sum(tril(p4<0.01,-1)))/tot_pairs;
fprintf('Mean fractions sig = %0.2f, %0.2f, %0.2f, %0.2f\n', frac1, frac2, frac3, frac4);

% correlation heatmaps
figure('Position', [200, 200, 1300, 250]);
subplot(141)
imagesc(r1.*(p1<0.01));caxis([-1 1]);title('WM Ref');
colormap(RB_spect_cmap); xticks(40:40:size(mat1, 2)); yticks(40:40:size(mat1, 2));
subplot(142)
imagesc(r2.*(p2<0.01));caxis([-1 1]);title('Neutral Ref');
colormap(RB_spect_cmap); xticks(40:40:size(mat1, 2)); yticks(40:40:size(mat1, 2)); yticklabels([]);
subplot(143)
imagesc(r3.*(p3<0.01));caxis([-1 1]);title('WM-Rereferenced');
colormap(RB_spect_cmap); xticks(40:40:size(mat1, 2)); yticks(40:40:size(mat1, 2)); yticklabels([]);
subplot(144)
imagesc(r4.*(p4<0.01));caxis([-1 1]);title('Neut-Rereferenced');
colormap(RB_spect_cmap); xticks(40:40:size(mat1, 2)); yticks(40:40:size(mat1, 2)); yticklabels([]);

saveas(gcf, fullfile('output', sprintf('hm_corr_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('hm_corr_%s', dataStr)), 'svg');

% Save a colorbar
figure('Position', [400, 400, 200, 100]);
imagesc(r4.*(p4<0.01));caxis([-1 1]);title('Neut-Rereferenced')
colormap(RB_spect_cmap)
colorbar;
saveas(gcf, fullfile('output', 'hm_corr_colorbar'), 'svg');

%% 3) Calculate latencies of positive, negative responses by peakfinder
% Code from EB, mostly unchanged

for x=1:ndatasets

    if x==1
        data = mat1;

        % prepare empty arrays for latency, amplitude, prominence (pos and neg peaks) at the first iteration
        P1_lat=zeros(size(data,2), ndatasets);
        P1_amp=zeros(size(data,2), ndatasets);
        P1_prom=zeros(size(data,2), ndatasets);
        N1_lat=zeros(size(data,2), ndatasets);
        N1_amp=zeros(size(data,2), ndatasets);
        N1_prom=zeros(size(data,2), ndatasets);

    elseif x==2
        data = mat2;
    elseif x==3
        data = mat3;
    elseif x==4
        data = mat4;
    end
    
    % search between 9 and 100 ms post stim: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10076215/
    time_start = find(time >= 9.0000 & time<9.5);
    time_end = find(time >=100 & time<100.5);
    data_cut = data(time_start:time_end, :);
    time_cut = time(time_start:time_end);
    
    % remove some high-freq noise to smooth out the EP trace
    fs=2000;
    y = highpass(data_cut,300,fs);
    data_cl = data_cut - y;
    
    % plot sanity check
    %nch=20;
    %figure;plot(time_cut,data_cut(:,nch))
    %hold on
    %plot(time_cut,data_cl(:,nch))
    
    % baseline (unused for now)
    base_startend = find(time >= -110.0000 & time<120.000);
    base_cut = data(base_startend, :);
        
    % search for N1 and P1 
    minprom = 20; % this could be adaptive based on pre-stim standard dev?
    % min(std(base_cut))/10 or something of the sort. the range of the
    % re-ref and neutral dataset is much lower than the original, would be
    % nice to have min prominence be based on that, but at the same time
    % for consistency we will not do that and just keep it the same across
    % datasets

    for nch=1:size(data,2)
    
        [maxPKS,maxLOCS,~,maxPROM] = findpeaks(data_cl(:,nch),"MinPeakProminence",minprom);
        [minPKS,minLOCS,~,minPROM] = findpeaks(-data_cl(:,nch),"MinPeakProminence",minprom);
        
        if ~isempty(maxLOCS)
            % if there are peaks detected, get the one with the highest
            % prominence
            if numel(maxPROM)>1
                clear sort_ind
                [~,sort_ind] = sort(maxPROM,'descend');
                P1_lat(nch,x) = time_cut(maxLOCS(sort_ind(1))); 
                P1_amp(nch,x) = maxPKS(sort_ind(1));
                P1_prom(nch,x) = maxPROM(sort_ind(1));
            else
                P1_lat(nch,x) = time_cut(maxLOCS); 
                P1_amp(nch,x) = maxPKS;
                P1_prom(nch,x) = maxPROM;
            end
        else
            P1_lat(nch,x) = NaN; 
            P1_amp(nch,x) = NaN;
            P1_prom(nch,x) = NaN;
        end
        
        if ~isempty(minLOCS)
            if numel(minPROM)>1
                clear sort_ind
                [~,sort_ind] = sort(minPROM,'descend');
                N1_lat(nch,x) = time_cut(minLOCS(sort_ind(1))); 
                N1_amp(nch,x) = -minPKS(sort_ind(1));
                N1_prom(nch,x) = -minPROM(sort_ind(1));
            else
                N1_lat(nch,x) = time_cut(minLOCS); 
                N1_amp(nch,x) = -minPKS;
                N1_prom(nch,x) = -minPROM;
            end
        else
            N1_lat(nch,x) = NaN; 
            N1_amp(nch,x) = NaN;
            N1_prom(nch,x) = NaN;
        end

    end

end

%% 4) plot histograms of N1, P1 latencies

% colors are edited separately
% https://colorbrewer2.org/#type=sequential&scheme=Blues&n=3
% Blues: from light to dark
Blues = [222 235 247;158 202 225;49 130 189;]./255;
Greens = [229, 245, 224; 161, 217, 155; 49, 163, 84]./255;
Reds = [254 224 210; 252 146 114; 222 45 38]./255;
Purples = [239, 237, 245; 188, 189, 220; 117, 107, 177]./255;
bin_width_ms = 5;

% Negative peaks (N1s)
figure;
histogram(N1_lat(:,1), 0:bin_width_ms:100,'FaceColor',Blues(1,:), Normalization="probability")
hold on
histogram(N1_lat(:,2), 0:bin_width_ms:100,'FaceColor',Greens(1,:), Normalization="probability")
histogram(N1_lat(:,3), 0:bin_width_ms:100,'FaceColor',Blues(3,:), Normalization="probability")
histogram(N1_lat(:,4), 0:bin_width_ms:100,'FaceColor',Greens(3,:), Normalization="probability")
xlim([0, 100]); ylim([0 0.8]);
saveas(gcf, fullfile('output', sprintf('N1_latencies_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('N1_latencies_%s', dataStr)), 'svg');

% Positive peaks (P1s)
figure;
histogram(P1_lat(:,1), 0:bin_width_ms:100,'FaceColor',Reds(1,:), Normalization="probability")
hold on
histogram(P1_lat(:,2), 0:bin_width_ms:100,'FaceColor',Purples(2,:), Normalization="probability")
histogram(P1_lat(:,3), 0:bin_width_ms:100,'FaceColor',Reds(3,:), Normalization="probability")
histogram(P1_lat(:,4), 0:bin_width_ms:100,'FaceColor',Purples(3,:), Normalization="probability")
xlim([0, 100]); ylim([0 0.8])
saveas(gcf, fullfile('output', sprintf('P1_latencies_%s', dataStr)), 'png');
saveas(gcf, fullfile('output', sprintf('P1_latencies_%s', dataStr)), 'svg');

%% 5) Calculate mean, mode frequency, and entropy of N1, P1 distributions

fprintf('Order: contaminated, neutral, contaminated CARLA, neutral CARLA\n');

% Mean +/- SD
means_N1 = mean(N1_lat, 'omitnan');
sd_N1 = std(N1_lat, 'omitnan');
means_P1 = mean(P1_lat, 'omitnan');
sd_P1 = std(P1_lat, 'omitnan');
fprintf('N1 latency means: '); fprintf('%f, ', means_N1); fprintf('\n');
fprintf('N1 latency SD: '); fprintf('%f, ', sd_N1); fprintf('\n');
fprintf('P1 latency means: '); fprintf('%f, ', means_P1); fprintf('\n');
fprintf('P1 latency SD: '); fprintf('%f, ', sd_P1); fprintf('\n');

% Mode, within 1 ms
N1_lat_round = round(N1_lat);
P1_lat_round = round(P1_lat);
modes_N1 = mode(round(N1_lat));
modes_P1 = mode(round(P1_lat));
modes_N1_pct = 100*sum(modes_N1 == N1_lat_round) ./ sum(~isnan(N1_lat_round));
modes_P1_pct = 100*sum(modes_P1 == P1_lat_round) ./ sum(~isnan(P1_lat_round));
fprintf('N1 mode size: '); fprintf('%f, ', modes_N1_pct); fprintf('\n');
fprintf('P1 mode size: '); fprintf('%f, ', modes_P1_pct); fprintf('\n');

% calculate entropy of distributions from binned probabilities
entropyfun = @(x) -sum(x.*log2(x));
nboot = 10000; % number of samples to bootstrap for SD
rng(0);

% get bootstrap indices (bootstrap channels)
[~, bootsam] = bootstrp(nboot, [], N1_lat(:, 1));

% calculate entropy for each bootstrapped sample
entropy_N1 = zeros(nboot, 4);
entropy_P1 = zeros(nboot, 4);

for bb = 1:nboot

    % 20 5ms bins
    binned_N1 = zeros(20, 4);
    binned_P1 = zeros(20, 4);

    % bootstrap channels
    N1_lat_boot = N1_lat(bootsam(:, bb), :);
    P1_lat_boot = P1_lat(bootsam(:, bb), :);

    for ii = 1:20
        for jj = 1:4
            binned_N1(ii, jj) = sum(N1_lat_boot(:, jj) >= (ii-1)*5 & N1_lat_boot(:, jj) < ii*5);
            binned_P1(ii, jj) = sum(P1_lat_boot(:, jj) >= (ii-1)*5 & P1_lat_boot(:, jj) < ii*5);
        end
    end

    % normalize bins by total number
    binned_N1_norm = binned_N1 ./ sum(binned_N1);
    binned_P1_norm = binned_P1 ./ sum(binned_P1);

    % store entropy value
    for ii = 1:4
        entropy_N1(bb, ii) = entropyfun(binned_N1_norm(binned_N1_norm(:, ii) > 0, ii));
        entropy_P1(bb, ii) = entropyfun(binned_P1_norm(binned_P1_norm(:, ii) > 0, ii));
    end
end

entropy_N1_mean = mean(entropy_N1);
entropy_P1_mean = mean(entropy_P1);
entropy_N1_SD = std(entropy_N1);
entropy_P1_SD = std(entropy_P1);
fprintf('N1 entropy means: '); fprintf('%f, ', entropy_N1_mean); fprintf('\n');
fprintf('N1 entropy SD: '); fprintf('%f, ', entropy_N1_SD); fprintf('\n');
fprintf('P1 entropy means: '); fprintf('%f, ', entropy_P1_mean); fprintf('\n');
fprintf('P1 entropy SD: '); fprintf('%f, ', entropy_P1_SD); fprintf('\n');

if strcmpi(dataStr, 'baylor'), return; end

%% 6) Calculate mean, mode frequency, and entropy of response durations (Mayo sub2 only)
% CRP response duration histogram already saved by comparePEP_SFN2024

load(fullfile('output', sprintf('CRP_durations_%s.mat', dataStr)));

CRP_durs = [dursW, dursG, dursWCarla, dursGCarla];

% Mean +/- SD
means_CRP = mean(CRP_durs, 'omitnan');
sd_CRP = std(CRP_durs, 'omitnan');
fprintf('Duration means: '); fprintf('%f, ', means_CRP); fprintf('\n');
fprintf('Duration SD: '); fprintf('%f, ', sd_CRP); fprintf('\n');

% Mode, rounded to nearest 10 ms
CRP_durs_round = round(CRP_durs, 2);
modes_CRP = mode(CRP_durs_round);
modes_CRP_pct = 100*sum(modes_CRP == CRP_durs_round) ./ sum(~isnan(CRP_durs));
fprintf('Duration mode size: '); fprintf('%f, ', modes_CRP_pct); fprintf('\n');

% calculate entropy of distributions from binned probabilities
entropyfun = @(x) -sum(x.*log2(x));
nboot = 10000; % number of samples to bootstrap for SD
rng(0);

% get bootstrap indices (bootstrap channels)
[~, bootsam] = bootstrp(nboot, [], CRP_durs(:, 1));

% calculate entropy for each bootstrapped sample
entropy_CRP = zeros(nboot, 4);

for bb = 1:nboot

    % 20 50ms bins
    binned_CRP = zeros(20, 4);

    % bootstrap channels
    CRP_durs_boot = CRP_durs(bootsam(:, bb), :);

    for ii = 1:20
        for jj = 1:4
            binned_CRP(ii, jj) = sum(CRP_durs_boot(:, jj) >= (ii-1)*0.05 & CRP_durs_boot(:, jj) < ii*0.05);
        end
    end

    % normalize bins by total number
    binned_CRP_norm = binned_CRP ./ sum(binned_CRP);

    % store entropy value
    for ii = 1:4
        entropy_CRP(bb, ii) = entropyfun(binned_CRP_norm(binned_CRP_norm(:, ii) > 0, ii));
    end
end

entropy_CRP_mean = mean(entropy_CRP);
entropy_CRP_SD = std(entropy_CRP);
fprintf('Duration entropy means: '); fprintf('%f, ', entropy_CRP_mean); fprintf('\n');
fprintf('Duration entropy SD: '); fprintf('%f, ', entropy_CRP_SD); fprintf('\n');
