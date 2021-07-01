function [fid_out,traj_out] = remove_extra_points(fid_in,traj_in,Method_Params)

%In this function, remove zero filling and acquisition shift points.

%% Remove zero filling:
numZeroFill = find(squeeze(fid_in(:,1,1,1,1,1))==0);
numZeroFill(numZeroFill <10)=[]; %Sometimes the first couple of points are 0, but we don't want to delete those
fid_in(numZeroFill,:,:,:,:,:)=[];

%% Remove AcqShift Points
traj_in(1:Method_Params.AcqShift,:,:,:,:,:) = [];
fid_in(:,1:Method_Params.AcqShift,:) = [];

fid_out = fid_in;
traj_out = traj_in;