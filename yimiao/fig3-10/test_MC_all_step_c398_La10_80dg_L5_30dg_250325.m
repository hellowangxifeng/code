%% MC与2Dbsf 理论值 在接收平面第一象限直线y=tan(phi)上对比
clear
close
set(0,'defaultfigurecolor','w');
Color = {[0,0,0],[1,0.5,0],[1,0.8,0.25],[0.53,0.84,0.20],[0.5,0.75,0.9],[1,0.6,0.6],[0.5,0.4,0.8]};
%% figure
% 理论
P0=1; La=10; div_rad=0.01;  N_times=1; water_type='c0.398'; theta0=80*pi/180;
phi0=30*pi/180;
% Lw=div_rad*La/cos(theta0+div_rad)*sin(phi0);
Lw=5;
r=linspace(0,0.5,100);
fi=[0 30 45 60 90]*pi/180;
bsf_2D=zeros(length(fi),length(r));
y_m1=zeros(length(fi),length(r));
r_m1=zeros(length(fi),length(r));
tic
for j =1:length(fi)
    x=r*cos(fi(j)); y=r*sin(fi(j));
    mu=asin(La*sin(pi/2-theta0)./(sqrt((y/cos(phi0)).^2+La^2-2*(y/cos(phi0))*La*cos(pi/2-theta0))));
    for i=1:length(r)
        j,i
        y_m1(j,i)=y(i)+(Lw+y(i)*tan(phi0))/tan(pi-theta0-mu(i));
        r_m1(j,i)=sqrt(y_m1(j,i)^2+x(i)^2);
        bsf_2D(j,i)= M1_bsf_240729_7_2(P0,La,Lw,theta0,phi0,div_rad,water_type,x(i),y(i));
    end
end
time=toc
save('org_RISM2_c398_La10_80dg_L5_30dg_240825_7_2.mat','r','r_m1','bsf_2D','fi','time','theta0','phi0');
% load RIS_c151_La75_80dg_M1_80dg_240704

%% MC
lambda=532e-9;
Np=5e7;
n_MC=10; 
fi=[0 30 45 60 90]*pi/180;
rx_0_end=[];rx_30_end=[];rx_45_end=[];rx_60_end=[];rx_90_end=[];
tic
for i =1:n_MC
i
N_times=1;
% [para_temp,aaa_M]=MC_all_step(P0,La,L0,div_rad,theta0,N_times,water_type);%发射平面
[para_temp,aaa] = MC_all_step_240809(P0,La,Lw,div_rad,theta0,phi0,N_times,Np,water_type);
% [MCpara] = MC_a_2_w_in_sea_aixs1_2D_UWOC_240623(P0, 1, Np,water_type, L, div_rad, lambda,pi/2-theta0,0*pi);
% load MC_rPlane_c151_La10_pi3_M1_30dg_240704
% load MC_bsf_c151_L15_pi4_240701
% load MC_bsf2D_c151_L10_pi4_240623
% load aaa
num_photon=Np;
aaa_end = aaa;
rx_0 =  getbsf_in_fi_240617(aaa_end,fi(1),num_photon);
rx_30 = getbsf_in_fi_240617(aaa_end,fi(2),num_photon);
rx_45 = getbsf_in_fi_240617(aaa_end,fi(3),num_photon);
rx_60 = getbsf_in_fi_240617(aaa_end,fi(4),num_photon);
rx_90 = getbsf_in_fi_240617(aaa_end,fi(5),num_photon);
rx_0_end =[rx_0_end;rx_0];
rx_30_end =[rx_30_end;rx_30];
rx_45_end =[rx_45_end;rx_45];
rx_60_end =[rx_60_end;rx_60];
rx_90_end =[rx_90_end;rx_90];
end
time=toc
rx_0_out  = mean(rx_0_end,1);
rx_30_out=mean(rx_30_end,1);
rx_45_out=mean(rx_45_end,1);
rx_60_out=mean(rx_60_end,1);
rx_90_out=mean(rx_90_end,1);
rx_out=[rx_0_out;rx_30_out;rx_45_out;rx_60_out;rx_90_out];
save('MC_bsf_fi_c398_L10_80dg_L5_30dg_240825.mat','r','rx_out','fi','time');
%%   c=0.151 Lw=15  
figure
load org_RISM2_c398_La10_80dg_L5_30dg_240825_7_2
load MC_bsf_fi_c398_L10_80dg_L5_30dg_240825
bsf=bsf_2D;  %%理论值
mc=rx_out;  %%MC  
r=linspace(0,0.5,100);
semilogy(-1, 10, 'x','Color',Color{1},'LineWidth', 1.2 );hold on
for i=1:5
    semilogy(r, bsf(i,:) , '-','Color',Color{1+i},'LineWidth', 1.8 );hold on
end
for i=1:5
    mc1=mc(i,:);
    semilogy(r(1:10:100), mc1(1:10:100), 'x','Color',Color{1},'LineWidth', 1.2 );hold on
end
line1 = '{\fontname{华文中宋}理论}  \phi=0^°';line2='{\fontname{华文中宋}仿真}';
line3 = '{\fontname{华文中宋}理论}  \phi=30^°';line4='{\fontname{华文中宋}仿真}  \phi=30^°';
line5 = '{\fontname{华文中宋}理论}  \phi=45^°';line6='{\fontname{华文中宋}仿真}  \phi=45^°';
line7 = '{\fontname{华文中宋}理论}  \phi=60^°';line8='{\fontname{华文中宋}仿真}  \phi=60^°';
line9 = '{\fontname{华文中宋}理论}  \phi=90^°';line10='{\fontname{华文中宋}仿真}  \phi=90^°';
Le=legend({line2,line1,line3,line5,line7,line9},'box','off','NumColumns', 2);
% legend({line1,line3,line5,line7,line9,line2},'box','off','NumColumns', 2);
% title('$c$=0.398, $L$=10m, $\theta=80^\circ$','interpreter','latex');
XL=xlabel('$r$(m) ','interpreter','latex');
YL=ylabel(' $I_{M}({r},{\phi}$)(W$\cdot$m$^{-2}$)','interpreter','latex');
ylim([2e-3 2e-1])
set(gca,'YTick',[1e-3,1e-2,1e-1,1e0,1e1]);
ax = gca;
set(ax,'FontSize',18,'FontName','Times New Roman');
set(Le,'FontSize',18,'FontName','Times New Roman');
set(XL,'FontSize',20,'FontName','Times New Roman');
set(YL,'FontSize',20,'FontName','Times New Roman');
pbaspect([1.268 1 1])