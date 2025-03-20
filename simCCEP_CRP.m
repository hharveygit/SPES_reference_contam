%% Script to simulate a CCEP/BSEP to label positive, negative peaks, and response duration in Figure 2
%
%    If this code is used in a publication, please cite the manuscript:
%    "Proper reference selection and re-referencing to mitigate bias in single pulse electrical stimulation data"
%    by H. Huang, J.A. Adkinson, M.A. Jensen, M. Hasen, I.A. Danstrom, K.R. Bijanki, N.M. Gregg, K.J. Miller,
%    S.A. Sheth, D. Hermes, E. Bartoli.
%
%    simCCEP_CRP.m
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
%% 

set(0, 'DefaultFigureRenderer', 'painters');

srate = 600;
tt = -0.5:1/srate:1-1/srate;

ntrs = 12; % number of trials

%% A) Create data, add stim artifacts and EP

% find a good model BSEP
seed = 27;
rng(seed);

% stores the data
V0 = zeros(length(tt), ntrs);

Aart = 20 + rand(ntrs, 1)*4;
artifact = sin(2*pi*200*tt)';
artifact(tt < 0 | tt > 0.005) = 0;
V0 = V0 + artifact*Aart';

A = 100;
sig = genRandSig(tt, 1, A);

% shift initial signal a bit in time
sig = circshift(sig, 0.01*srate);

V1 = V0 + sig;

figure('Position', [200, 200, 400, 200]); plot(tt, V1);
xlim([-0.1, 0.3]);
xlabel('Time (s)'); ylabel('Channels');

%% B) add random brown noise across all trials

rng(2*seed);

noiseRand = cumsum(randn(3*length(tt), ntrs), 1); % give triple the number of time points so we can highpass it
noiseRand = ieeg_highpass(noiseRand, srate, true);
noiseRand = noiseRand((length(tt)+1) : 2*length(tt), :);
V2 = V1 + noiseRand;

figure('Position', [200, 200, 400, 200]); hold on
plot(tt, V2, 'Color', [0.5, 0.5, 0.5]);
plot(tt, mean(V2, 2), 'k-', 'LineWidth', 1.5);
xlim([-0.1, 0.5]);
xlabel('Time (s)'); ylabel('Channels');

%% Find negative and postive peaks using same way as real data

% positive peak
figure; findpeaks(mean(V2(tt >= 0.009 & tt < 0.1, :), 2), "MinPeakProminence", 20, 'Annotate', 'extents');

% negative peak
figure; findpeaks(-mean(V2(tt >= 0.009 & tt < 0.1, :), 2), "MinPeakProminence", 20, 'Annotate', 'extents');

%% Determine CRP duration
[crp_parms, crp_projs] = CRP_method(V2(tt >= 0.009, :), tt(tt >= 0.009));

figure('Position', [200, 200, 400, 200]); hold on
plot(tt, V2, 'Color', [0.5, 0.5, 0.5]);
plot(tt, mean(V2, 2), 'k-', 'LineWidth', 1.5);
xline(crp_parms.tR, 'Color', 'r');
xlim([-0.1, 0.5]);
saveas(gcf, fullfile('output', 'metrics_demo'), 'png');
saveas(gcf, fullfile('output', 'metrics_demo'), 'svg');
