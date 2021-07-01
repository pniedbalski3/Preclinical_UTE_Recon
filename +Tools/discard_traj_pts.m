function [FID_out,traj_out] = discard_traj_pts(FID_mat,traj)

%Function to shape trajectories for Jim Pipe's Recon:
%Actually, already in the correct shape, just need to make sure values are
%between +/- 0.5

rad = squeeze(sqrt(traj(1,:,:).^2+traj(2,:,:).^2+traj(3,:,:).^2));

toobig = rad>0.5;

%I need to keep everything the right shape, so find the lowest point that
%needs to be discarded.
pts = size(rad,1)+1;
for i = 1:size(traj,3)
    toobig = nnz(rad(:,i)>0.5);
    if toobig < pts 
        pts = toobig;
    end
end

if pts < size(rad,1)+1
    FID_mat((end-pts+1):end,:,:,:,:,:) = [];
    traj(:,(end-pts+1):end,:) = [];
end

traj_out = traj;
FID_out = Fid_mat;
    


