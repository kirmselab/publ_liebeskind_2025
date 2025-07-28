%% This script computes changes in Ihold and Rm in voltage-clamp recordings.
% ====================================================================================================
% Authors: Lissy Liebeskind, Knut Kirmse
% last modification: 10/2023
% ====================================================================================================
% How to use:
% (1) Specify parameters
% (2) Run the script
% (3) De-select bad traces in app
% ====================================================================================================
% Inputs:
% (1) TXT file (episodic stimulation)
        % col1: time (ms)
        % col2: current (pA); col3: voltage (mV) > sweep1
        % col4: current (pA); col5: voltage (mV) > sweep2 etc.
% ====================================================================================================
% Outputs:
% (1) MAT file (in the same dir as the TXT file)
% ====================================================================================================
%
%% Specify parameters
opt.path = 'C:\MyFolder\'; % path
opt.file = 'testdata.txt'; % TXT file name
opt.fs = 50000; % sampling frequency [Hz]
%
% --- define Ihold periods ---
% test pulse is from 1.0 s to 1.1 s; we allow Ihold to settle until 1.2 s
opt.a = 0;      % start of first Ihold period (before test pulse) [s]
opt.b = 1;      % end of first Ihold period (before test pulse) [s]
opt.c = 1.2;    % start of second Ihold period (after test pulse) [s]
opt.d = 10;     % end of second Ihold period (after test pulse) [s]
%
% --- define test pulse (periods) for Rm measurements ---
opt.testpulse = -10;            % test pulse amplitude [mV]
opt.hold_range_on  = 0.958;     % onset of time range to compute Ihold [s]; default: 0.990;
opt.hold_range_off = 0.998;     % offset of time range to compute Ihold [s]
opt.peak_range_on  = 0.998;     % onset of time range to compute Ipeak [s]
opt.peak_range_off = 1.010;     % offset of time range to compute Ipeak [s]
opt.ss_range_on    = 1.058;     % onset of time range to compute Iss [s]; default: 1.090;
opt.ss_range_off   = 1.098;     % offset of time range to compute Iss [s]
%
%% Clean-up
clearvars -except opt; close all;
opt.timestamp = datestr(now, 'yymmdd_HHMMSS');
%
%% Compute Ihold per sweep
% convert time from seconds to sample points
a = opt.a * opt.fs + 1;         
b = opt.b * opt.fs;
c = opt.c * opt.fs + 1;
d = opt.d * opt.fs;
hold_range = int32((opt.hold_range_on * opt.fs + 1):(opt.hold_range_off * opt.fs));
peak_range = int32((opt.peak_range_on * opt.fs + 1):(opt.peak_range_off * opt.fs));
ss_range   = int32((opt.ss_range_on * opt.fs + 1):(opt.ss_range_off * opt.fs));
%  read data
disp('Loading file ...')
data = readmatrix([opt.path, opt.file]);
disp(['File loaded: ', opt.path, opt.file])
Nsweeps = (size(data, 2)-1)/2;
if data(2, 1) ~= 1/opt.fs*1000
    error('Sampling interval is incorrect.')
end
% compute Ihold
data_Ihold = data([a:b, c:d], 2:2:end); % data used for estimating Ihold (i.e. excl. test pulse)
%
Hist_range = [floor(min(data_Ihold(:))), ceil(max(data_Ihold(:)))];
Hist_edges = Hist_range(1):0.5:Hist_range(2);
Hist_bincenters = movmean(Hist_edges, [0, 1], 'Endpoints', 'discard');
%
Median_Ihold_perSweep = median(data_Ihold);
Hist_Ihold_perSweep = NaN(numel(Hist_edges)-1, Nsweeps);
for k = 1:Nsweeps
    Hist_Ihold_perSweep(:, k) = histcounts(data_Ihold(:, k), Hist_edges, 'Normalization', 'probability');
    if abs(sum(Hist_Ihold_perSweep(:, k))-1) > 1e-10
        disp(['For trace ', num2str(k), ', the sum of relative frequencies per bin does not equal zero.'])
    end
end
[~, Mode_Ihold_perSweep] = max(Hist_Ihold_perSweep, [], 1);
Mode_Ihold_perSweep = Hist_bincenters(Mode_Ihold_perSweep);
clear k
% compute Rm and Rs
% calculate absolute values of I_hold, I_ss and I_peak 
I_ss =  mean(data(ss_range,2:2:end)); 
I_hold =  mean(data(hold_range,2:2:end));
if opt.testpulse > 0
    I_peak = max(data(peak_range,2:2:end));
else
    I_peak = min(data(peak_range,2:2:end));
end
% calculate Ra Rm and baseline
deltaIpeak = I_peak - I_hold;
deltaIss = I_ss - I_hold;
Ra = opt.testpulse ./ deltaIpeak;
Rm = (opt.testpulse ./ deltaIss ) - Ra;
Ra = Ra .* 1000; % [MOhm]
Rm = Rm .* 1000; % [MOhm]
%
%% Plot in ViewerApp
appexit = 0;
TonicCurrent_Viewer
while appexit == 0
    pause(1)
end
%
%% Complete & save MAT file
VarNames = {'Sweep number []', 'Modal Ihold (pA)', 'Selected traces Ihold []', 'Rs [MOhm]', 'Rm [MOhm]', 'Selected traces Rm []'};
SweepNumber = (1:Nsweeps)';
Results = table(SweepNumber, Mode_Ihold_perSweep', selectedTraces_Ihold', Ra', Rm', selectedTraces_Rm');
Results.Properties.VariableNames(:) = VarNames;
clear appexit data data_Ihold SweepNumber VarNames
%
save([opt.path, opt.file, '_Results', '_', opt.timestamp, '.mat'])
disp(['File saved: ', opt.path, opt.file, '_Results', '_', opt.timestamp, '.mat'])
disp('Completed.')
%