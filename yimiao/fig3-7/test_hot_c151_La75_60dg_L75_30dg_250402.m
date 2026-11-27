%% MC与2Dbsf 理论值 在接收平面第一象限直线y=tan(phi)上对比
clear
close
set(0,'defaultfigurecolor','w');
Color = {[0,0,0],[1,0.5,0],[1,0.8,0.25],[0.53,0.84,0.20],[0.5,0.75,0.9],[1,0.6,0.6],[0.5,0.4,0.8]};
Colors = [  2  48  71
      18 104 131
      39 158 188
     144 201 230
     255 255 255]/255;
nColors = 256;
newmap = interp1(1:size(Colors,1), Colors, linspace(1, size(Colors,1), nColors), 'linear');
%% figure
% MC
P0=1;  div_rad=0.01; water_type='c0.151'; theta0=60*pi/180;
La=7.5; Lw=7.5;
phi0=30*pi/180; phi1=0*pi/180;
lambda=532e-9;
Np=5e7; 
 
tic
N_times=1;
% [para_temp,aaa_M]=MC_all_step(P0,La,L0,div_rad,theta0,N_times,water_type);%发射平面
dmax=0.3; m1_g=200;
x_dir = linspace(-dmax, dmax,m1_g+1) ;
y_dir = linspace(-dmax, dmax,m1_g+1) ;

[para_temp,aaa] = MC_all_step_x_250319_5(P0,La,Lw,div_rad,theta0,phi0,phi1,N_times,Np,water_type);

num_photon=Np;
aaa_end = aaa;
h_out = hotdotpic_240618(x_dir,y_dir,aaa_end.x',aaa_end.y',aaa_end.w,Np,dmax,m1_g); 
time=toc

save('hot_bsf_c151_L75_60dg_L75_30_250402.mat','x_dir','y_dir','h_out','dmax','m1_g','Np','time');
%%
figure
load hot_bsf_c151_L75_60dg_L75_30_250402
colormap(newmap);  
imagesc(x_dir(1:end-1)+dmax/m1_g,y_dir(1:end-1)+dmax/m1_g,h_out); 
XL=xlabel('$y_M$(m)','interpreter','latex');
YL=ylabel('$x_M$(m)','interpreter','latex');
ylim([-dmax dmax]);
xlim([-dmax dmax]);
% title('$c$=0.151, $L$=10m, $\theta=60$','interpreter','latex');
colorbar;
ax=gca;
set(ax,'XTick',[(-0.3:0.1:0.3)]);%x轴刻
set(ax,'FontSize',20,'FontName','Times New Roman');
% set(Le,'FontSize',18,'FontName','Times New Roman');
set(XL,'FontSize',22,'FontName','Times New Roman');
set(YL,'FontSize',22,'FontName','Times New Roman');
% pbaspect([1.268 1 1])

