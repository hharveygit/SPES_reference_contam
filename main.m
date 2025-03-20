%% Code walkthrough to generate all results and figures
%
%    If this code is used in a publication, please cite the manuscript:
%    "Proper reference selection and re-referencing to mitigate bias in single pulse electrical stimulation data"
%    by H. Huang, J.A. Adkinson, M.A. Jensen, M. Hasen, I.A. Danstrom, K.R. Bijanki, N.M. Gregg, K.J. Miller,
%    S.A. Sheth, D. Hermes, E. Bartoli.
%
%    Main script
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
% Required packages:
% - mnl_ieegBasics: https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
% - CRP_method is already included in ./functions to avoid separate installation. See file for original reference.
%
% Download the data files from OSF (https://osf.io/ryfqt/) and copy them to a "data" folder in the working directory (./data)

% set paths, make output directory
addpath('functions');
addpath('data');
mkdir('output');
addpath('output');

%% Run modified CARLA on Baylor and Mayo datasets
% Saves all figure 1 panels/results,
% as well as intermediates, CAR waveforms, and response duration histogram (sub 2) for figure 2
% Note that CARLA is not deterministic and outputs will vary slightly on repeat iteration

dataStr = 'baylor'; % to configure which subject is asved
comparePEP_SFN2024

dataStr = 'mayo_LK1-LK2';
comparePEP_SFN2024

%% Rest of results for figure 2

% This saves correlation heatmaps, positive/negative peak latency histograms, values in Table S1
dataStr = 'baylor';
Latency_FindPeaks_EB

dataStr = 'mayo_LK1-LK2';
Latency_FindPeaks_EB

% This saves the model PEP with response duration marked in Figure 2, bottom left
simCCEP_CRP
