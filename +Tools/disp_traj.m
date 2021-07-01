function disp_traj(traj,slow)

if nargin < 2
    slow = false;
end

figure('Name','Trajectory Sanity Check')
view(-41.1,16.2);
hold on
for i = 1:size(traj,3)
    plot3(traj(1,:,i),traj(2,:,i),traj(3,:,i))
    if slow
        pause(0.1)
    end
end
hold off

