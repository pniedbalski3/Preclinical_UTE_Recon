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
%Tools.disp_traj(traj(:,:,:,1),0);


%% Now, Reconstruct Data
%We can do either Pipe or Robertson Recon Methods: Either way, everything
%is reconstructed with the same trajectories, so calculate the density
%compensation separately before looping through and regridding:
%% Density Compensation
disp('Begin Density Compensation')
if contains(Method_Params.Dims,'3D')
    nIter = 10; %10 generally seems a good number of iterations for density compensation
    DCF = Recon.get_DCF(traj,Method_Params.MatrixSize,nIter);
    Image = zeros([Method_Params.MatrixSize size(FID_mat,3) size(FID_mat,4) size(FID_mat,5) size(FID_mat,6) size(FID_mat,7)]);
    Fin_Image = zeros([Method_Params.MatrixSize size(FID_mat,3) size(FID_mat,4) size(FID_mat,5) size(FID_mat,7)]);
else
    nIter = 10;
    traj = Tools.column_traj(traj);
    DCF = Tools.get_DCF_Robertson(traj,Method_Params.MatrixSize,nIter);
    Image = zeros([Method_Params.MatrixSize(1) Method_Params.MatrixSize(2) size(FID_mat,3) size(FID_mat,4) size(FID_mat,5) size(FID_mat,6) size(FID_mat,7)]);
    Fin_Image = zeros([Method_Params.MatrixSize(1) Method_Params.MatrixSize(2) size(FID_mat,3) size(FID_mat,4) size(FID_mat,5) size(FID_mat,7)]);
end
disp('Finished Density Compensation')
%% Reconstruction
disp('Begin Reconstruction')
for i = 1:size(FID_mat,3)
    for j = 1:size(FID_mat,4)
        for k = 1:size(FID_mat,5)
            for m = 1:size(FID_mat,7)
                for l = 1:size(FID_mat,6)
                    fid_tmp = squeeze(FID_mat(:,:,i,j,k,l,m));
                   % figure;imagesc(abs(fid_tmp))
                    if contains(Method_Params.Dims,'3D')
                        Image(:,:,:,i,j,k,l,m) = Recon.pipe_recon(Method_Params.MatrixSize(1),fid_tmp,traj,DCF,l,size(FID_mat,6));
                    else
                        fid_recon = reshape(fid_tmp,1,[])';
                        Image(:,:,i,j,k,l,m) = DCF.reconstruct(fid_recon, traj);
                    end
                end
                if contains(Method_Params.Dims,'3D')
                    Fin_Image(:,:,:,i,j,k,m) = Tools.soscoilcombine(squeeze(Image(:,:,:,i,j,k,:,m)));
                else
                    Fin_Image(:,:,i,j,k,m) = Tools.soscoilcombine(squeeze(Image(:,:,i,j,k,:,m)));
                end
            end
        end
    end
end
Fin_Image = squeeze(Fin_Image);
disp('Finish Reconstruction')
                



