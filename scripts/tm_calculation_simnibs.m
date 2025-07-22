% Assuming hm_840_simnibs and elec_840 have already been created as per your script
% and loaded into the workspace

% Define the subset of channels you want to include
channels = {'Fp1', 'Fpz', 'Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2',...
'FC6','M1','T7','C3','Cz','C4','T8','M2','CP5','CP1','CP2','CP6','P7','P3',...
'Pz','P4','P8','POz','O1','Oz','O2','AF7','AF3','AF4','AF8','F5','F1','F2',...
'F6','FC3','FCz','FC4','C5','C1','C2','C6','CP3','CPz','CP4','P5','P1','P2',...
'P6','PO5','PO3','PO4','PO6','FT7','FT8','TP7','TP8','PO7','PO8'};

% EMG1 EMG2 and Trigger are not in this 10-10 EEG system. So we have 64
% channels, instead of the 67 channels that we had previously in the EEG
% analysis.

elec_selected = select_elec_channels(elec, channels);

% Create a new configuration structure for ft_prepare_headmodel
cfg = [];
cfg.method = 'simbio';
cfg.conductivity = [0.126 0.275 1.654 1e-6 0.465 0.5 0.008 0.025 0.6 0.16];
cfg.tissuelabel = {'white', 'gray', 'csf', 'holes', 'scalp', 'eyeballs', 'compactbone', 'spongybone', 'blood', 'muscle'}'; 


cfg.elec = elec_selected;
cfg.headmodel = headmodel_simnibs;

% Use ft_prepare_headmodel with the specified channels
headmodel_selected = ft_prepare_headmodel(cfg);

% Plot the head model and the selected electrodes
figure;
ft_plot_headmodel(headmodel_selected, 'vertexcolor', 'skin_medium');
hold on;
ft_plot_sens(elec_selected, 'elecshape', 'sphere', 'label', 'on');
camlight;

%% Compute the transformation matrix
[tm_simnibs_defaced, sens_simnibs_defaced] = ft_prepare_vol_sens(headmodel_selected, elec_selected);

% Save the selected electrode configuration and head model
% save tm_840_simnibs tm_840_simnibs -v7.3
% save sens_840_simnibs sens_840_simnibs -v7.3

% Function to select specified channels from the electrode structure
function elec_selected = select_elec_channels(elec, selected_channels)
    % Find the indices of the selected channels
    [sel_idx, ~] = find(ismember(elec.label, selected_channels));
    
    % Ensure indices are valid
    sel_idx = sel_idx(sel_idx > 0);
    
    % Create a new electrode structure with only the selected channels
    elec_selected = [];
    elec_selected.label = elec.label(sel_idx);
    elec_selected.chanpos = elec.chanpos(sel_idx, :);
    elec_selected.elecpos = elec.elecpos(sel_idx, :);
    
    if isfield(elec, 'unit')
        elec_selected.unit = elec.unit;
    end
end