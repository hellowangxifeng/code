% set(0,'defaultfigurecolor','w');
% P0=1;
% La=10;
% Lw=5;
% div_rad=0.01;
% theta0=60/180*pi;
% fi=30/180*pi;
% N_times=1;
% water_type='c0.151';
% Np=2e7;
% Lw=div_rad*La/cos(theta0)*sin(fi);
function [para_temp,aaa] = MC_all_step_x_250319_5(P0,La,Lw,div_rad,theta0,fi,fi0,N_times,Np,water_type)
%% z_air2sea_allstep 带有大气MC
%% Copy_of_z_air2sea_allstep 没有大气MC
% P0=1;
% div_rad = 0.001; La=100;
w0=div_rad*La;
% theta_all=[pi/2,pi/4,pi/3,pi/6];%与y轴夹角

%MC 空中
% visibility = 10; 
lambda =532e-9;
% N_times=10;
aaa=[];

%     theta=theta_all(i);
    theta = theta0;
    La_real = La;
 for j =1:N_times
%                 j
        [para_1] = MC_a_2_w_in_sea_aixs1_2D_UWOC_240710_2(P0, 1,Np, water_type, La, div_rad, lambda,theta0,0);
%         [para_1] = MC_a_2_w_in_sea_aixs1(P0, 1, visibility, La_real, div_rad, lambda,theta);  %   已和初始的大气MC对比 没有什么问题  出来的是海平面坐标系下坐标
        % [para_1] = no_aMC_getinfo(theta,div_rad,La);%注意此处的权重w设置的都为1
        % [rec_x,rec_y,rec_w] = z_MC_a_2_w_all(P0, 1, visibility, 100, div_rad, lambda);
        para_temp = para_1;
%         sum(para_temp.w)/1e7
aaa_M=para_1;
% dmax=0.2; m1_g=200;
% x_dir = linspace(-dmax, dmax,m1_g+1) ;
% y_dir = linspace(-dmax, dmax,m1_g+1) ;
% h_out = hotdotpic_240618(x_dir,y_dir,aaa_M.x',aaa_M.y',aaa_M.w,Np,dmax,m1_g); 
% save('MC_ris_c151_La10_pi3_240704.mat','x_dir','y_dir','para_1','h_out','dmax','m1_g','Np');
% load newmap
% % load MC_ris_c151_L5_pi3_240703
% figure
% colormap(newmap);  
% imagesc(x_dir(1:end-1)+dmax/m1_g,y_dir(1:end-1)+dmax/m1_g,h_out); 
% xlabel('{\ity} [m]');
% ylabel('{\itx} [m]');
% ylim([-0.15 0.15]);
% xlim([-0.15 0.15]);
% title('$c$=0.151, $La$=10m, $\theta=60^\circ$','interpreter','latex');
% colorbar;
% ax=gca;
% set(ax,'FontSize',14,'FontName','Times New Roman');
%%
%折射
% V10 = 5;
% n_sea=get_n_only_1(V10,1); %抽样得到  一定在YOZ平面
% clear
% aaa=[];
% P0=1;
% La=10;
% theta0=60/180*pi;
% fi=30/180*pi;
% div_rad=0.01;
% 
% N_times=1;
% water_type='c0.151';
% Np=2e7;
% load MC_ris_c151_L5_pi3_240703
        para_1.uz=-para_1.uz;
        ux_new1=para_1.ux';
        uy_new1=para_1.uy'*cos(fi-theta0)-para_1.uz'*sin(fi-theta0);
        uz_new1=para_1.uy'*sin(fi-theta0)+para_1.uz'*cos(fi-theta0);
        ux_new2=ux_new1*cos(fi0)+uz_new1*sin(fi0);
        uy_new2=uy_new1;
        uz_new2=-ux_new1*sin(fi0)+uz_new1*cos(fi0);
%         n_sea=[0,0,-1];
%         [ux_new1 ,uy_new1,uz_new1 ] = getnewdir_afterRIS_240704(para_1.ux',para_1.uy',para_1.uz',n_sea,theta0,fi); %折射  新的方向向量已经进行归一化
        para_1.ux=ux_new1; para_1.uy=uy_new1; para_1.uz=uz_new1; %
        para_1.ux=ux_new2; para_1.uy=uy_new2; para_1.uz=uz_new2;
        % para_1.ux = (para_1.ux)'; para_1.uy = (para_1.uy)';  para_1.uz = (para_1.uz)';    %无折射时要添加这一行 有折射时不用添加
        % sum(para_1.w)/1e7
        
%由海平面坐标系转到水下MC坐标系
%         sin_fi = 0.75*sin(pi/2 - theta); %入射角为pi/2 - theta
%         fi=asin(sin_fi);
        j_3 =[0,cos(fi),-sin(fi)];             k_3 = [0,sin(fi),cos(fi)];         i_3 = cross(j_3,k_3);
        j_3 =j_3/norm(j_3);                 k_3 =k_3/norm(j_3);             i_3 =i_3/norm(j_3);  
        T13=[i_3;j_3;k_3];
        k_4 = [-sin(fi0),0,cos(fi0)];         i_4 =[cos(fi0),0,sin(fi0)];      j_4 = [0,1,0];   
     j_4 =j_4/norm(k_4);                 k_4 =k_4/norm(k_4);             i_4 =i_4/norm(k_4);
     T14=[i_4;j_4;k_4];
        uxyz_3 = [(para_1.ux)', (para_1.uy)', (para_1.uz)']*T14*T13;    %方向向量转换   
        xyz_3   = [(para_1.x), (para_1.y), (para_1.z)]*T14*T13;          %位置坐标转换  之前没转换 最后结果不对 看看是不是这块的原因

%         uxyz_3 = [(para_1.ux)', (para_1.uy)', (para_1.uz)']*T13;    %方向向量转换   
%         xyz_3   = [(para_1.x), (para_1.y), (para_1.z)]*T13;          %位置坐标转换  之前没转换 最后结果不对 看看是不是这块的原因

        para_3.ux = (uxyz_3(:,1))';  para_3.uy =(uxyz_3(:,2))' ;  para_3.uz = (uxyz_3(:,3))'; 
        para_3.x  =  xyz_3(:,1);       para_3.y  = xyz_3(:,2);       para_3.z  = xyz_3(:,3);   %xyz uxuyuz都为列向量
        para_3.w  = para_1.w;

%水下MC仿真
%         water_type='c0.151';
%计算发射平面M的光强信息
%         d_left = La*sin(div_rad)/sin(theta-div_rad);
%         sin_fi = 0.75*sin(pi/2 - theta); %入射角为pi/2 - theta
        % L0=d_left*sin_fi;   %%L0太小了 在MC仿真中起不到作用？
%         Lw=15;

        [aaa1] = MC_sea_in_3axis_250319_2(P0, 1, water_type,Lw ,para_3.x,para_3.y,para_3.w,para_3.ux,para_3.uy,para_3.uz,fi0);   %189注释
        aaa=[aaa;aaa1];
        
% dmax=0.2; m1_g=200;
% x_dir = linspace(-dmax, dmax,m1_g+1) ;
% y_dir = linspace(-dmax, dmax,m1_g+1) ;
% h_out = hotdotpic_240618(x_dir,y_dir,aaa1.x',aaa1.y',aaa1.w,Np,dmax,m1_g); 
% save('MC_rPlane_c151_La10_pi3_M1_30dg_240704.mat','x_dir','y_dir','aaa1','h_out','dmax','m1_g','Np','aaa1');
% load MC_rPlane_c151_La10_pi3_M1_30dg_240704
% load newmap
% figure
% colormap(newmap);  
% imagesc(x_dir(1:end-1)+dmax/m1_g,y_dir(1:end-1)+dmax/m1_g,h_out); 
% xlabel('{\ity} [m]');
% ylabel('{\itx} [m]');
% ylim([-0.2 0.2]);
% xlim([-0.2 0.2]);
% title('$c$=0.151, $L_a$=10m, M1, $\theta_0=60^\circ$, $\phi=30^\circ$','interpreter','latex');
% colorbar;
% ax=gca;
% set(ax,'FontSize',14,'FontName','Times New Roman');
%%
end
    
end

