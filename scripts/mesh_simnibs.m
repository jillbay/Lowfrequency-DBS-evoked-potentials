%% Add simnibs -> matlab path
addpath('C:\Users\Jill\SimNIBS-4.1\simnibs_env\Lib\site-packages\simnibs\matlab_tools'); 

addpath('C:\Users\Jill\Documents\UTwente\Master\fieldtrip-20240326');
ft_defaults;
%% Load mesh
mesh = mesh_load_gmsh4('C:\Users\Jill\Documents\UTwente\Master\Preliminary results\DBS014\Intermediate Results\m2m_DBS014\DBS014.msh');
ft_plot_mesh(mesh)
%% Rewrite to fieldtrip format
simnibs_mesh = struct(); 

simnibs_mesh.pos = mesh.nodes; 
simnibs_mesh.tet = mesh.tetrahedra; 
simnibs_mesh.tissue = mesh.tetrahedron_regions; 
simnibs_mesh.unit = 'mm'; 
simnibs_mesh.coordsys = 'ras';
simnibs_mesh.tissuelabel = {'white', 'gray', 'csf', 'holes', 'scalp', 'eyeballs', 'compactbone', 'spongybone', 'blood', 'muscle'}'; 

%% Reorientate tetrahedral 
% The following error is displayed if you try to create the headmodel like
% this:
% error('Elements have wrong orientation, consider exchanging node 3 and 4');

% So you have to exchange node 3 and 4: 
% this function reorders nodes in a mesh to ensure correct orientation
simnibs_mesh.tet = meshreorient(simnibs_mesh.pos, simnibs_mesh.tet);
%%
figure; 
ft_plot_mesh(simnibs_mesh)
%% Create headmodel
cfg = [];
cfg.method = 'simbio'; 
cfg.conductivity = [0.126 0.275 1.654 1e-6 0.465 0.5 0.08 0.025 0.6 0.16]; % from Simnibs
cfg.tissuelabel = {'white', 'gray', 'csf', 'holes', 'scalp', 'eyeballs', 'compactbone', 'spongybone', 'blood', 'muscle'}';
cfg.coordsys = 'ras'; 
headmodel_simnibs = ft_prepare_headmodel(cfg, simnibs_mesh);
% the call to "ft_prepare_headmodel" took 796 seconds and required the additional allocation of an estimated 929 MB
 %%
figure;
ft_plot_headmodel(headmodel_simnibs, 'vertexcolor', 'skin_medium');

%% Load and Align Electrodes
elec_online = ft_read_sens('standard_1020.elc');
cfg = [];
cfg.method = 'interactive';
cfg.elec = elec_online;
cfg.headshape = headmodel_simnibs;
elec = ft_electroderealign(cfg);

% DBS011
% [0 0 0]
% [ 0.99 1.02 1]
% [-1 86 14]

% DBS009
% [5 -10 0]
% [1 0.96 0.95]
% [0 38 -3]

% DBS014
% [5 -6 0]
% [0.95 1 1]
% [-1 85 20]

%% Save and Plot
% elec_840_simnibs = elec;
figure;
ft_plot_headmodel(headmodel_simnibs, 'vertexcolor', 'skin_medium');
hold on;
ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'on');
camlight;
%%
elec_842_simnibs = elec_840_simnibs; 
hm_842_simnibs = hm_840_simnibs;

%%
save headmodel_simnibs headmodel_simnibs -v7.3
save elec elec -v7.3


