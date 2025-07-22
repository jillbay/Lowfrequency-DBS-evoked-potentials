%% Artifact removal with template subtraction and linear interpolation as done by Peeters et al. 
% This script will remove the artifact from the low-frequency DBS datasets
% s000820 of patient 7

%% Load Fieldtrip
addpath('C:\Users\Jill\Documents\UTwente\Master\fieldtrip-20240326');
ft_defaults;
%% Load the DBS file
addpath('C:\Users\Jill\Documents\UTwente\Master\EANSKE\DBS007\EEG');
filename = 's0000820a'; % define the data path and its name

% Read events:
cfg            = [];
cfg.dataset    = [filename, '.edf'];            % set the name of the dataset
cfg.continuous = 'yes';                         % continuous signal
cfg.channel    = 'all';                         % define channel type
data           = ft_preprocessing(cfg);         % read raw data
%% Rereference to common average
cfg = [];
cfg.reref = 'yes';
cfg.refmethod = 'avg';
cfg.refchannel = 'all';
cfg.demean        = 'yes';
data = ft_preprocessing(cfg, data);
%% Check spectogram (Optional)
% figure                                                                      % FIGURE 1, spectogram
% pspectrum(data.trial{1}(19,:), data.fsample, 'spectrogram', 'FrequencyLimits', [10 600], 'TimeResolution',1);   % Creates powerspectrum of channel 19 in this case
% title(['Dataset-', filename])
%% Create powerspectrum (Optional)
% window = ceil(5 * data.fsample); %5 seconds window
% nfft = 2^(2+nextpow2(window)); 
% noverlap=0.5*window; %define amount of overlap between windows
% fr = data.fsample/nfft; %Frequency resolution
% [pxx_on,f_on] = pwelch(data.trial{1}(19,:),window,noverlap,nfft,data.fsample); %Calculate power spectrum
% figure, plot(f_on,pxx_on) % Plot power spectrum
%% Apply montages
myelec_label= {'Fp1', 'Fpz', 'Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2',...
'FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3',...
'Pz','P4','P8','POz','O1','Oz','O2','AF7','AF3','AF4','AF8','F5','F1','F2',...
'F6','FC3','FCz','FC4','C5','C1','C2','C6','CP3','CPz','CP4','P5','P1','P2',...
'P6','PO5','PO3','PO4','PO6','FT7','FT8','TP7','TP8','PO7','PO8','EMG1','EMG2','Trigger'}'; % Names of channel labels  (electrodes)

montage = [];                                  % label montage
montage.tra      = eye(size(data.label,1));
montage.labelold = data.label;
montage.labelnew = myelec_label;

data    = ft_apply_montage(data,montage);
%% Check channels visually
cfg          = [];
cfg.viewmode = 'vertical'; %Channels below each other
cfg.channel  = {'M2','M1','FC3'}; %Select channels you want to check; can scroll afterwards through others
 ft_databrowser(cfg, data);   
%% use M2 to define trials
chan4peaks = data.trial{1}(strcmp(data.label, 'M2'),:); %Select channel M2
figure, plot(data.time{1},chan4peaks) %Plot signal of M2
xlabel('time (s)')
ylabel('mV')
title('M2')
[peaks, locs] = findpeaks(chan4peaks, 'MinPeakDistance', 25, 'MinPeakHeight', 1); %Select peaks based on M2
figure, plot(chan4peaks), hold on, plot(locs, peaks, 'r+')  
%% Select all DBS peaks
[peaks_all, locs_all] = findpeaks(chan4peaks,'MinPeakDistance',5,'MinPeakHeight',2);
figure, plot(chan4peaks), hold on, plot(locs_all, peaks_all, 'r+') %20Kpeaks
%% keep only the train with lower frequency
locs_middle = locs_all > 8296760-1 & locs_all < 8534745+1;
locs = locs_all(locs_middle);
peaks = peaks_all(locs_middle);

% seperate odd and even peaks as stimulation in left and right hemisphere
locs_odd = locs(1:2:end); 
peaks_odd = peaks(1:2:end);

locs_even = locs(2:2:end);
peaks_even = peaks(2:2:end);

figure, plot(chan4peaks), hold on, plot(locs_odd, peaks_odd, 'r+'), hold on, plot(locs_even, peaks_even, 'g+') %20Kpeaks
%% keep only the train with high frequency (optional, keep another train)
% locs_middle = locs_all > 8591460 & locs_all < 9018420;
% locs = locs_all(locs_middle);
% peaks = peaks_all(locs_middle);
% 
% locs_odd = locs(1:2:end);
% peaks_odd = peaks(1:2:end);
% 
% locs_even = locs(2:2:end);
% peaks_even = peaks(2:2:end);
% 
% figure, plot(chan4peaks), hold on, plot(locs_odd, peaks_odd, 'r+'), hold on, plot(locs_even, peaks_even, 'g+') %20Kpeaks
%% keep only the first train (Optional to keep another train)
% locs_middle = locs_all > 270915 & locs_all < 741311;
% locs = locs_all(locs_middle);
% peaks = peaks_all(locs_middle);
% 
% locs_odd = locs(1:2:end);
% peaks_odd = peaks(1:2:end);
% 
% locs_even = locs(2:2:end);
% peaks_even = peaks(2:2:end);
% 
% figure, plot(chan4peaks), hold on, plot(locs_odd, peaks_odd, 'r+'), hold on, plot(locs_even, peaks_even, 'g+') %20Kpeaks
%% reject "artifacts"
cfg        = [];
cfg.metric = 'range'; 
cfg.method = 'summary'; % use by default summary method
data       = ft_rejectvisual(cfg,data); %  POz, P2, TP7, PO8, EMG1, EMG2, M1
%% decide pre and post
msBefore = 30; %Define time before peak
msAfter = 250-31; %Define time after peak
pretrl = msBefore*4; 
posttrl = msAfter*4; 

trl_startstop = [];

%Select epoch times of the odd peaks
trl_startstop(:,1) = locs_odd-pretrl; 
trl_startstop(:,2) = locs_odd+posttrl; 

%Select epoch times of the even peaks
% trl_startstop(:,1) = locs_even-pretrl;
% trl_startstop(:,2) = locs_even+posttrl; 

cfg_trl = [];
cfg_trl.trl = [trl_startstop -pretrl*ones(size(trl_startstop,1),1)];
data = ft_redefinetrial(cfg_trl, data);

trl = cfg_trl.trl; 

%% timelock analysis
cfg = [];
cfg.preproc.demean = 'yes';
cfg.preproc.baselinewindow = [-0.030 -0.005]; %Low freq
cfg.covariance       = 'yes';
avgERPdbs64_raw = ft_timelockanalysis(cfg, data);
%% reject "artifacts"
cfg        = [];
cfg.metric = 'range';  
cfg.method = 'summary'; % use by default summary method
avgERPdbs64_raw       = ft_rejectvisual(cfg,avgERPdbs64_raw); %  Trigger
%% Check epoched data
figure, plot(1000*avgERPdbs64_raw.time, 1000*avgERPdbs64_raw.avg(:,:)) %Plot all channels; Big artifact visible
xlabel('time(ms)'), ylabel('Amplitude (uV)')
title('DBS artifact of the odd peaks')
ylim([-10 10]); 
xlim([-30 220]);
grid on; 
%%
% figure, plot(1000*avg_820.time, 1000*avg_820.avg(:,:)) %Plot all channels; Big artifact visible
% xlabel('time(ms)'), ylabel('Amplitude (uV)')
% title('DBS artifact of the odd peaks')
% ylim([-10 10]); 
% xlim([-30 220]);
% grid on; 
%% Save dataset
avg_820=avgERPdbs64_raw;
% save('avg_820_odd_nofiltering', 'avg_820');
% save('avg_820_even_nofiltering', 'avg_820');
%% Template subtraction based on mastoid electrodes
sig = avgERPdbs64_raw.avg;  % Select the epoched signal
tim = avgERPdbs64_raw.time; % Select timevector

% Averaging M1 (13) and M2 (19) (template)
M2     = sig(18,123:end);      % M1 is removed
template  = M2;                 

%% Calculate the relative error
RE_all = zeros(length(sig),length(template));  % Predefine relative error matrixrt

for i = 1:length(sig(:,1))
    for j = 1:length(template)
        RE_all(i,j) = ((template(1,j)-sig(i, 123+j-1))) ./ template(1,j); %Calculation relative errors; the RE is calculated for every timestep in the signal
         if isnan(RE_all(i,j))  %If template and sig are exactly similar; put RE on 1
            RE_all(i,j)= 1;
         end
    end
    avgRE(i)    = mean(RE_all(i,2:5),2); %Calculate average over RE 2:5; we have a scaling factor for every channel left
end
%% CHECK THE TEMPLATE
% Check template with regards to mastoid electrodes
figure, plot(M2), 
hold on, plot(template) %Check fitting of M1, M2 and template in 1 figure
xlabel('Time(ms)'), ylabel('Amplitude(uV)') 
legend('M1','M2','Template')
title('Created template')

% Check the template with regards to the whole signal
figure
plot(1000*avgERPdbs64_raw.time,1000*sig,'color',[0,0,0.5]); %Plot whole epoched signal
hold on, 
plot(1000*avgERPdbs64_raw.time(123:end),1000*template,'color',[1,0,0],'linewidth',3); %Plot template from n=123
 xlabel('Time(ms)'), ylabel('Amplitude(uV)') 
%% Fit & subtract the template
fitted_signal = zeros(length(sig(:,1)),length(template));  %Predefine matrix; This is amount of channels x length of template
fitted_templ  = zeros(length(sig(:,1)),length(template));  %Predefine matrix; We create a template for every channel

chan=length(avgERPdbs64_raw.label); %Amount of channels
for i = 1:chan 
    for j = 1:length(template) %Go through every timestep in the template
        fitted_templ(i,j)    = (1-avgRE(1,i)) * template(1,j); %Rescale every channel with its own individual template based on 1-RE
        fitted_signal(i,j)   = sig(i, 123+j-1) - fitted_templ(i,j); %Subtract template from signal
    end
end

figure, plot(1000*tim(123:end),1000*sig(59, 123:end),'LineWidth',3), hold on %Plot EEG signal
plot(1000*tim(123:end),1000*fitted_templ(59, :),'LineWidth',2.5), hold off  %Plot fitted template to compare the fit
title('Fitting of template with FP1')
xlabel('Time(ms)'), ylabel('Amplitude(uV)')  
ylim([-20 20])

figure, plot(1000*tim(123:end), 1000*fitted_signal(59,:),'LineWidth',2),
xlabel('Time(ms)'), ylabel('Amplitude(uV)')  
title('Signal after artifact removal with mastoid templates')
xlabel('Time(ms)'), ylabel('Amplitude(uV)') 

%% Apply linear interpolation until n=123
% Create matrix with second part of signal and untouched first part of signal (A)
A            = zeros(size(sig)); %Predefine matrix
A(:,1:123)    = sig(:,1:123); %Untouched first part of signal, we will apply linear interpolation on this
A(:,123:end)  = fitted_signal(:,1:end);  %Signal after template subtraction. n=14 for better results
 
% The goal here is to perform linear interpolation on the initial segment
% of each channel (n=1:123), based on the signal's mean amplitude after
% artifact removal, effectively smoothing the transition between the raw
% and corrected segments. 

% Interpolation for all epochs: starting in 0
interp_i = [];
for i = 1:chan
    x_i  = [tim(1), tim(123)]; %Select x axis
    A_mean(i)=mean(A(i,123:end),2); %Calculate the mean of the 2nd part of the signal. If it is the first value it interpolates to an outlier
    v_i  = [0, A_mean(i)];        % Interpolation between y value 0 and the previously calculated mean
    % defines the y-coordinates for interpolation, starting from 0 to the
    % calculated mean amplitude. 
    xq_i = tim(1):2.5000e-04:tim(123); %Steps in x direction
    % defines the query points for interpolation along the x-axis, covering
    % the time span from 0:123 with steps determined by the sampling
    % interval 2.5e-4 seconds. 
    vq_i     = interp1(x_i, v_i, xq_i); % interpolation 
    interp_i = [interp_i; vq_i]; % interpolated segments matrix. 
end

%% Add interpolation to filtered signal
% A            = zeros(size(sig));
% A(:, 123:end) = fitted_signal(:,1:end);
 A(:,1:123)    = interp_i;  %Add interpolation to first part of the signal

 
 figure, plot(1000*tim,1000*A(:,:))

 xlabel('Time(ms)'), ylabel('Amplitude(uV)')
 title('Signal after artifact removal & interpolation')
 legend()
 ylim([-20 20])

 avgERPdbs64_raw.avg = A; %Rename to its old name
%% Save filtered data
avg_820=avgERPdbs64_raw;
%save('avg_820_odd_filtered','avg_820')
%save('avg_820_even_filtered','avg_820')
%% Global Mean Field
cfg = [];
cfg.method = 'amplitude';
EEG_free_gmfp = ft_globalmeanfield(cfg, avgERPdbs64_raw); %Calculate absolute average
figure, 
plot(1000*avgERPdbs64_raw.time,1000*avgERPdbs64_raw.avg,'color',[0,0,0.5]); %All channels
hold on;
plot(1000*avgERPdbs64_raw.time,1000*EEG_free_gmfp.avg,'color',[1,0,0],'linewidth',3);  %Check signal with global mean field. The blue lines are all channels, 
%The red line is the absolute average of all channels
xlabel('time(ms)'), ylabel('uV')
title('Global mean field after artifact removal')
ylim([-10  10])

avg_odd=avgERPdbs64_raw;
%save('avg_820')
 %% Topoplot 
% addpath('Z:\General\LindaLoosveld\Final\Intermediate results')
% load layout_m10_v3.mat
cfg = [];
cfg.layout = layout_m10; %Electrode locations
layout = ft_prepare_layout(cfg);
ft_layoutplot(cfg);

% Prompt the user to enter a time point in milliseconds
prompt = 'Enter a time point (in milliseconds) for the topoplot: ';
selected_time_ms =6; %Time of topoplot (ms)
selected_time = selected_time_ms / 1000;  % Convert from milliseconds to seconds

% Find the closest time point index to the selected time
[~, idx] = min(abs(avgERPdbs64_raw.time - selected_time));
toi = avgERPdbs64_raw.time(idx);

[mxx, idxm] = max(max(abs(avgERPdbs64_raw.avg(:, idx))));
toi_mean_trial = toi(idxm);

cfg = [];
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.xlim = [toi_mean_trial, toi_mean_trial + 0.01 * toi_mean_trial];
cfg.layout = layout;
% cfg.marker = 'labels'; %Add electrode names in plot
cfg.fontsize = 14;
cfg.zlim = [-5 5]; %Defines color axis
figure;
ft_topoplotER(cfg, avgERPdbs64_raw);
colorbar;

time_title=1000*toi_mean_trial;
title(['Topoplot at time ', num2str(time_title), ' ms']);
%% Topoplot video
load layout_m10_v3.mat
cfg = [];
cfg.layout = layout_m10';
layout = ft_prepare_layout(cfg);
ft_layoutplot(cfg);

% Define the time points for the topoplots
times2check = linspace(0, max(avgERPdbs64_raw.time)); % Times to plot

% Create a VideoWriter object
videoFile = 'topoplot_video.avi';
v = VideoWriter(videoFile, 'Motion JPEG AVI');
open(v);

% Interpolate between the time points for smoother transition
interpolated_times = linspace(min(avgERPdbs64_raw.time), max(avgERPdbs64_raw.time), 50); %Adjust 50 to amount of frames in the video

for t = 1:length(interpolated_times)
        [~, idx] = min(abs(avgERPdbs64_raw.time - interpolated_times(t))); % Find the closest time point index to the desired time point
    toi = avgERPdbs64_raw.time(idx);
    
    [mxx, idxm] = max(max(abs(avgERPdbs64_raw.avg(:, idx))));
    toi_mean_trial = toi(idxm);
    
    cfg = [];
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.xlim = [toi_mean_trial, toi_mean_trial + 0.01];
    cfg.layout = layout;
    cfg.fontsize = 14;
   % cfg.marker='labels'   
    figure;
    ft_topoplotER(cfg, avgERPdbs64_raw); %Show frames when generating
    colorbar;
    title(['Topoplot at time ', num2str(toi_mean_trial), ' seconds']);
    
    % Capture the current figure as a frame for the video
    frame = getframe(gcf);
    writeVideo(v, frame.cdata);
    
    % Pause to display the frame for a longer duration
    pause(1); % Adjust the duration (in seconds) as needed
    
    close(gcf); % Close the figure to prevent accumulation of open figures
end

% Close the video writer
close(v);

% Display a message with the saved video file path
disp(['Video saved to: ', videoFile]);