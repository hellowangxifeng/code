

% this function postprocess wave gauge data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear

%% Para
% read each line of the text file
fid = fopen('Sine Freq=0.4Hz Amp=0.03m ang=0rad_20241103T110044.txt','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

%% Process the Wave Gauge
% eliminates the text in the file 
Data = A(1,6:end);
% loop the rest of the data
for i = 1:size(Data,2) - 1
    Data_char = Data{1,i};
    % extract the time stamp from the first few chars
    t_min = str2num(Data_char(4:5));
    t_sec = str2num(Data_char(7:12));
    t_cur = t_min*60 + t_sec;
    t_vec(i) = t_cur;
    % extact gauge data (six in current code)
    eta_gauge_temp = str2num(Data_char(13:end));
    eta_gauge_vec(i,:) = eta_gauge_temp(4:4:end);
end

%% Process Qualisys
load('FOWTmotion_11_3_30mmp4Hz.mat');
Qualisys_Data = Measurement28;
FOWT_3D_trans = squeeze(Qualisys_Data.RigidBodies.Positions);  
FOWT_3D_rot = squeeze(Qualisys_Data.RigidBodies.RPYs);  
dt = 1/Qualisys_Data.FrameRate;
t_vec_temp = 0:dt:(Qualisys_Data.Frames - 1)*dt;

%% Global time stamp calculation
% Calculate the global time stamp for Qualisys;  NOTE: Manual laps are
% utilized
lap2 = 2.01;  % 3.69
lap3 = 2.82;    % 5
t_vec_quali = t_vec_temp + lap3;
t_vec_gauge = t_vec - t_vec(1) + lap2;

%% Plot results
figure
plot(t_vec_gauge,eta_gauge_vec(:,1),'linewidth',1.5)
hold on
plot(t_vec_gauge,eta_gauge_vec(:,2),'linewidth',1.5)
plot(t_vec_gauge,eta_gauge_vec(:,3),'linewidth',1.5)
xlabel('Time')
ylabel('eta (m)')
legend('Gauge 1','Gauge 2','Gauge 3')
set(gca,'Fontsize',18,'FontWeight','bold')
grid on

figure
plot(t_vec_gauge,eta_gauge_vec(:,4),'linewidth',1.5)
hold on
plot(t_vec_gauge,eta_gauge_vec(:,5),'linewidth',1.5)
plot(t_vec_gauge,eta_gauge_vec(:,6),'linewidth',1.5)
xlabel('Time')
ylabel('eta (m)')
legend('Gauge 13','Gauge 7','Gauge 8')
set(gca,'Fontsize',18,'FontWeight','bold')
grid on

figure
plot(t_vec_gauge,eta_gauge_vec(:,6),'linewidth',1.5)
hold on
plot(t_vec_quali,(FOWT_3D_trans(3,:) - FOWT_3D_trans(3,1))/1000,'linewidth',1.5)
xlabel('Time')
ylabel('eta (m)')
legend('Gauge 8','FOWT heave')
set(gca,'Fontsize',18,'FontWeight','bold')
grid on

figure
plot(t_vec_quali,(FOWT_3D_rot(1,:) - FOWT_3D_rot(1,1)),'linewidth',1.5)
xlabel('Time')
ylabel('Rot (deg)')
set(gca,'Fontsize',18,'FontWeight','bold')
grid on