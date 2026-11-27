% Seabotix Stationkeeping Code
% 05/16/2016

% This code loads the analyzed data set for station keeping at
%   10 meters (P2) and at roughly 40 meters (P3). After loading the data
%   both the mean squared error (RMSE) and mean error (ME) are calculated
%   and a 3D plot of the track is generated.

% Dataset:
%   In each dataset are the follwing vectors:
%       -x: relative x position
%       -y: relative y position
%       -z: relative z position
%       -heading: heading of the vehicle
%       -t_stamp: time at which data recieved

%

%% P2 Part 1 - Station Keeping at 10 meters, run 1

clear
load('P2_stationkeeping_part1_final.mat');

% Find the length of all vectors
num_points = length(x);

% Vector to hold the error magnitude
err_mag = zeros(num_points, 1);

% Vector to hold the normalized heading (accounting for 0-360 wrap)
head_norm = zeros(num_points, 1);

% Loop through to calculate error / normalized heading
for i = 1:num_points
    err_mag(i) = sqrt(x(i)^2+y(i)^2+z(i)^2);
    head_norm(i) = min([abs(360 - heading(i)), abs(0 - heading(i))]);
end

% Generating the target - for this data it is trying to keep at [0, 0, 0]
target = zeros(num_points, 1);

% Find the RMSE and ME for overall, x, y, z and heading
[overall_RMSE, overall_me] = err_calc(err_mag,target);
[x_RMSE, x_me] = err_calc(x,target);
[y_RMSE, y_me] = err_calc(y,target);
[z_RMSE, z_me] = err_calc(z,target);
[heading_RMSE, heading_me] = err_calc(head_norm, target);

% Find the duration of station keeping
duration = t_stamp(end) - t_stamp(1);

% 3D plot of the Results
figure()
plot3(x,y,z,'b')
hold on
plot3(x(1),y(1),z(1),'g-x','LineWidth',3)
plot3(x(end),y(end),z(end),'r-o','LineWidth',3)
grid on
title('Station Keeping Results for P2 part 1')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')

% 3D compariosn plot of station keeping vs drifting
figure()
plot3(x,y,z,'b')
hold on
grid on
title('Station Keeping vs Drifting for P2 part 1')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')
load('seastate_P2.mat')
plot3(x,y,z,'r')
legend('Station Keeping','Drifting','Location','northeast')
view(-56,44)

% Bar graph of RMSE and ME
figure()

x_bars = [x_RMSE x_me];
y_bars = [y_RMSE y_me];
z_bars = [z_RMSE z_me];
overall_bars = [overall_RMSE overall_me];
heading_bars_1 = [heading_RMSE heading_me];

bar([x_bars; y_bars; z_bars; overall_bars])
set(gca,'XTickLabel',{'x error','y error','z error','overall error'})
legend('RMSE','ME','Location','northwest')
ylabel('Error (m)')
title('Station Keeping Results for P2 part 1')

figure()
bar(heading_bars_1)
set(gca,'XTickLabel',{'RMSE','ME'})
ylabel('Error (deg)')
title('Heading Error for P2 Part 1')

%% P2 Part 2 - Station Keeping at 10 meters, run 2

clearvars -except heading_bars_1 heading_bars_2 heading_bars_3
load('P2_stationkeeping_part2_final.mat');

% Find the length of all vectors
num_points = length(x);

% Vector to hold the error magnitude
err_mag = zeros(num_points, 1);

% Vector to hold the normalized heading (accounting for 0-360 wrap)
head_norm = zeros(num_points, 1);

% Loop through to calculate error / normalized heading
for i = 1:num_points
    err_mag(i) = sqrt(x(i)^2+y(i)^2+z(i)^2);
    head_norm(i) = min([abs(360 - heading(i)), abs(0 - heading(i))]);
end

% Generating the target - for this data it is trying to keep at [0, 0, 0]
target = zeros(num_points, 1);

% Find the RMSE and ME for overall, x, y, z and heading
[overall_RMSE, overall_me] = err_calc(err_mag,target);
[x_RMSE, x_me] = err_calc(x,target);
[y_RMSE, y_me] = err_calc(y,target);
[z_RMSE, z_me] = err_calc(z,target);
[heading_RMSE, heading_me] = err_calc(head_norm, target);

% Find the duration of station keeping
duration = t_stamp(end) - t_stamp(1);

% 3D plot of the Results
figure()
plot3(x,y,z,'b')
hold on
plot3(x(1),y(1),z(1),'g-x','LineWidth',3)
plot3(x(end),y(end),z(end),'r-o','LineWidth',3)
grid on
title('Station Keeping Results for P2 part 2')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')

% 3D compariosn plot of station keeping vs drifting
figure()
plot3(x,y,z,'b')
hold on
grid on
title('Station Keeping vs Drifting for P2 part 2')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')
load('seastate_P2.mat')
plot3(x,y,z,'r')
legend('Station Keeping','Drifting','Location','northeast')
view(-56,44)

% Bar graph of RMSE and ME
figure()

x_bars = [x_RMSE x_me];
y_bars = [y_RMSE y_me];
z_bars = [z_RMSE z_me];
overall_bars = [overall_RMSE overall_me];
heading_bars_2 = [heading_RMSE heading_me];

bar([x_bars; y_bars; z_bars; overall_bars])
set(gca,'XTickLabel',{'x error','y error','z error','overall error'})
legend('RMSE','ME','Location','northwest')
ylabel('Error (m)')
title('Station Keeping Results for P2 part 2')

figure()
bar(heading_bars_2)
set(gca,'XTickLabel',{'RMSE','ME'})
ylabel('Error (deg)')
title('Heading Error for P2 Part 2')

%% P3 Part 1 - Station Keeping at ~35 meters, run 1

clearvars -except heading_bars_1 heading_bars_2 heading_bars_3
load('P3_stationkeeping_part1_final.mat');

% Find the length of all vectors
num_points = length(x);

% Vector to hold the error magnitude
err_mag = zeros(num_points, 1);

% Vector to hold the normalized heading (accounting for 0-360 wrap)
head_norm = zeros(num_points, 1);

% Loop through to calculate error / normalized heading
for i = 1:num_points
    err_mag(i) = sqrt(x(i)^2+y(i)^2+z(i)^2);
    head_norm(i) = min([abs(360 - heading(i)), abs(0 - heading(i))]);
end

% Generating the target - for this data it is trying to keep at [0, 0, 0]
target = zeros(num_points, 1);

% Find the RMSE and ME for overall, x, y, z and heading
[overall_RMSE, overall_me] = err_calc(err_mag,target);
[x_RMSE, x_me] = err_calc(x,target);
[y_RMSE, y_me] = err_calc(y,target);
[z_RMSE, z_me] = err_calc(z,target);
[heading_RMSE, heading_me] = err_calc(head_norm, target);

% Find the duration of station keeping
duration = t_stamp(end) - t_stamp(1);

% 3D plot of the Results
figure()
plot3(x,y,z,'b')
hold on
plot3(x(1),y(1),z(1),'g-x','LineWidth',3)
plot3(x(end),y(end),z(end),'r-o','LineWidth',3)
grid on
title('Station Keeping Results for P3 part 1')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')

% Bar graph of RMSE and ME
figure()

x_bars = [x_RMSE x_me];
y_bars = [y_RMSE y_me];
z_bars = [z_RMSE z_me];
overall_bars = [overall_RMSE overall_me];
heading_bars_3 = [heading_RMSE heading_me];

bar([x_bars; y_bars; z_bars; overall_bars])
set(gca,'XTickLabel',{'x error','y error','z error','overall error'})
legend('RMSE','ME','Location','northwest')
ylabel('Error (m)')
title('Station Keeping Results for P3 part 1')

figure()
bar(heading_bars_3)
set(gca,'XTickLabel',{'RMSE','ME'})
ylabel('Error (deg)')
title('Heading Error for P3 Part 1')

%% P3 Part 2 - Station Keeping at ~35 meters, run 2

clearvars -except heading_bars_1 heading_bars_2 heading_bars_3
load('P3_stationkeeping_part2_final.mat');

% Find the length of all vectors
num_points = length(x);

% Vector to hold the error magnitude
err_mag = zeros(num_points, 1);

% Vector to hold the normalized heading (accounting for 0-360 wrap)
head_norm = zeros(num_points, 1);

% Loop through to calculate error / normalized heading
for i = 1:num_points
    err_mag(i) = sqrt(x(i)^2+y(i)^2+z(i)^2);
    head_norm(i) = min([abs(360 - heading(i)), abs(0 - heading(i))]);
end

% Generating the target - for this data it is trying to keep at [0, 0, 0]
target = zeros(num_points, 1);

% Find the RMSE and ME for overall, x, y, z and heading
[overall_RMSE, overall_me] = err_calc(err_mag,target);
[x_RMSE, x_me] = err_calc(x,target);
[y_RMSE, y_me] = err_calc(y,target);
[z_RMSE, z_me] = err_calc(z,target);
[heading_RMSE, heading_me] = err_calc(head_norm, target);

% Find the duration of station keeping
duration = t_stamp(end) - t_stamp(1);

% 3D plot of the Results
figure()
plot3(x,y,z,'b')
hold on
plot3(x(1),y(1),z(1),'g-x','LineWidth',3)
plot3(x(end),y(end),z(end),'r-o','LineWidth',3)
grid on
title('Station Keeping Results for P3 part 2')
xlabel('x distance (m)')
ylabel('y distance (m)')
zlabel('z distance (m)')

% Bar graph of RMSE and ME
figure()

x_bars = [x_RMSE x_me];
y_bars = [y_RMSE y_me];
z_bars = [z_RMSE z_me];
overall_bars = [overall_RMSE overall_me];
heading_bars_4 = [heading_RMSE heading_me];

bar([x_bars; y_bars; z_bars; overall_bars])
set(gca,'XTickLabel',{'x error','y error','z error','overall error'})
legend('RMSE','ME','Location','northwest')
ylabel('Error (m)')
title('Station Keeping Results for P3 part 2')

figure()
bar(heading_bars_4)
set(gca,'XTickLabel',{'RMSE','ME'})
ylabel('Error (deg)')
title('Heading Error for P3 Part 2')
