function preclinical_recon(path)

if nargin == 0
    path=uigetdir('C:\','Select Folder in Which Data is Stored');
end

%% Read in data
disp('Begin Reading Method File')
[traj,Method_Params] = Data_Import.read_method(path);
disp('Finished Reading Method File')

disp('Begin Reading Raw Data')
FID_mat = Data_Import.read_shape_bruker_data(path,Method_Params);
disp('Finished Reading Raw Data');

%Remove points that are known to be extraneous
[FID_mat,traj] = Tools.remove_extra_points(FID_mat,traj,Method_Params);
%Remove trajectory points that will misbehave
[FID_mat,traj] = Tools.discard_traj_pts(FID_mat,traj);

%% Optional Sanity Check
%Tools.disp_traj(traj,0);

%% Now, Reconstruct Data
%We can do either Pipe or Robertson Recon Methods: Either way, everything
%is reconstructed with the same trajectories, so calculate the density
%compensation separately before looping through and regridding:




