%% Source reconstruction on processed EEG data containing evoked potentials or epileptiform discharges 
% This script applies source reconstruction on dataset s000840 of patient DBS011.
%% Load Fieldtrip
addpath('C:\Users\Jill\Documents\UTwente\Master\fieldtrip-20221223\fieldtrip-20240326');
ft_defaults;
%% Load files to plot dipoles in MRI
addpath 'C:\Users\Jill\Documents\UTwente\Master\Preliminary results\DBS011\Intermediate Results'
%%
load leadfield_defaced.mat
load sens_simnibs.mat

%% Load MRI (If you want to plot dipoles to an MRI)
% addpath('Z:/General/LindaLoosveld/Final/Intermediate results')
 mri_orig = ft_read_mri('t1_anon.nii');
cfg = [];
ft_sourceplot(cfg,mri_orig);%%

%% Define specific timepoints in milliseconds
% dip = [8.25, 19.5, 48.5, 71.25]; %      left high C0
% dip = [9.5, 20, 58.5, 69.5] %         left low C0
% dip = [9, 20, 50.5, 68]; %             left high C1
% dip = [8.75, 20, 50, 67.75] %         left low C1
% dip = [7.25, 20, 40.25, 87.5] %       left high C2
% dip = [9.75, 22, 47.75, 72] %         left low C2
% dip = [8.75, 20.25, 35.5, 65.75] %    left high C3
% dip = [7.5, 16.75, 46.75, 77.5] %  left low C3

dip = 0; 

% Dipole fitting for specific timepoints
dip1 = cell(1, length(dip)); 

for ii = 1:length(dip)
    cfg_dipfit = [];
    cfg_dipfit.numdipoles = 1;
    cfg_dipfit.headmodel = []; % Must be empty as specified
    cfg_dipfit.sourcemodel = leadfield_defaced;
    cfg_dipfit.gridsearch = 'yes';
    cfg_dipfit.nonlinear = 'no';
    cfg_dipfit.elec = sens_simnibs;
    cfg_dipfit.latency = dip(ii) / 1000; 
    cfg_dipfit.channel = 'eeg';
    cfg_dipfit.unit = 'mm';
    cfg_dipfit.reducerank = 'no';

    source_eeg = ft_dipolefitting(cfg_dipfit, avg_820);
    dip1{1,ii} = source_eeg; 
end
%%
for i = 1:length(dip1)
    dipoles_840_left_low{4,i} = dip1{i}; 
end

% save dipoles_840_left_low dipoles_840_left_low -v7.3

%% Define base colors for each evoked potential
base_colors = [
    1 0 0;  % Red for EP1
    0 1 0;  % Green for EP2
    0 0 1;  % Blue for EP3
    1 1 0;  % Yellow for EP4
    1 0 1;  % magenta for EP5
    0 1 1;  % cyan for EP 6
];

figure;

hold on;

% Loop through each evoked potential
for ii = 1:length(dip1)
    num_dipoles = size(dip1, 1); 
    
    % Get the base color for this evoked potential
    base_color = base_colors(ii, :);
    
    % Generate shades of the base color for later EP dipoles
    for jj = 1:num_dipoles
        if ~isempty(dip1{jj, ii})

            shade_factor = 1 - (jj - 1) * 0.2; % Adjust the shade factor as needed
            color = base_color * shade_factor + [1 1 1] * (1 - shade_factor);
            
            % Plot the dipole
            ft_plot_dipole(dip1{jj, ii}.dip.pos(1,:), dip1{jj, ii}.dip.mom(1:3,:), 'color', color, 'unit', 'mm', 'alpha', 0.5);
        end
    end
end

% Plot MRI slices at the origin (if you want to add this)
hold on;
pos = [0, 45, 5]; % MRI crosses at (0,0,0)
ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1);
ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1);
ft_plot_slice(mri_orig.anatomy, 'transform', mri_orig.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1);

axis tight;
axis off;
hold off;


