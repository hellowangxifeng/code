% 2013_Monte-Carlo-Based Channel Characterization for Underwater Optical Communication Systems
% function [MC_BSF, Vz, RX_POW] = Photon_Propagation_Parallel(P0,number, water_type,psf_type,L, Rx_r, dr, div_rad, lambda)
%that is right
function [para_1,rx_pow_out,R0] = MC_a_2_w_in_sea_aixs1_2D_UWOC_240710_2(P0, number,Np, water_type, La, div_rad, lambda,theta0,theta02)
%输出已经是在海平面下的坐标系
%% generate CDF of phase scattering function\\分布函数
% [cdf_out, rad]= psf_cdf(psf_type);
% loglog(rad,cdf_out,'k-','LineWidth',2)   % 对比 simulation, modeling, and design of underwater optical communication systems, Fig. 4.27

% [a,b,c] = a_getabc(visibility); %得到大气中abc参数
[a,b,c] = Table_abc(water_type);%得到水下abc参数
% max_time = 0;
% min_time = 0;
speed = 3*10^8;  % speed of light in water %m/s
% R=Ra;
Timereach=[];
Wreach=[];
Wsum=[];
Wrow=cell(number,1);
%% 生成转换矩阵
% the
% k_1=[0,-cos(theta0),-sin(theta0)]; k_1 = k_1/norm(k_1);
% j_1=[0,1/cos(theta0)-cos(theta0),-sin(theta0)]; j_1 = j_1/norm(j_1); 
% i_1=cross(j_1,k_1); i_1=i_1/norm(i_1);
% T01 = [i_1;j_1;k_1]; 
% %%海平面方向向量[0，0，1]转换到发射机坐标系里
% n_0 = [0,0,1]*(T01^-1);
% a0=n_0(1); b0=n_0(2);c0=n_0(3); d0=-c0*La;
%% 生成转换矩阵结束
% Wsample=zeros(number,sample_point);
% rx_pow = zeros(number,length(Rx_r));
% mc_bsf = zeros(number, length(Rx_r));
    for nn=1:number
        %% Initialization
%         nn
%         Np = 1;
%         Np = 2*10^7;  % number of photons
        FOV = 180;
    %     div_deg = 1;  %0.01;  % divergence angle in degree
    %    div_rad = div_deg/180*pi; % divergence angle in radians
    %     lambda = 532e-9; % wavelength
    
        PreLoc = zeros(Np,3);  % pre-position
        PostLoc = zeros(Np,3);   % initial positions for (x,y,z)
        Wpre = P0*ones(1,Np);  % pre weight, Initially, each photon is launched into the medium with unit weight. 
        Wpost = P0*ones(1,Np);  % post weight,  Initially, each photon is launched into the medium with unit weight.
        time = zeros(1,Np);
        Dist = zeros(1,Np);    % Photon's distance traveled
        theta =zeros(1,Np); % polar angle in radians
        cosphi = zeros(1,Np);  % azimuthal angle in radians
        sinphi = zeros(1,Np);  % azimuthal angle in radians
        phi = zeros(1,Np);  % azimuthal angle in radians
        ScatNum  = zeros(1,Np);  % count scatter number for each photon
        
        PhStr = struct('Wpre',Wpre,'Wpost',Wpost, 'PreLoc', PreLoc, 'PostLoc',PostLoc,...
            'Dist', Dist, 'theta',theta,'phi',phi,'cosphi',cosphi,'sinphi',sinphi, 'ScatNum', ScatNum);  % photon characteristics
        
        Th = 1e-6;   % weight dropping threshold
        ah = 10; % rouletting threshold
        count = 1;
        curr_set = [1:Np];
            
        % Initiate scattering angles
        w0 =  lambda/(pi*div_rad);  % refer to :  https://en.wikipedia.org/wiki/Gaussian_beam#Radius_of_curvature,    Beam divergence
    
        rand_seed;
        wx=w0;
        wy=w0/cos(theta02);
        % phi = 2*pi*rand_init(2,:);
        % w0_phi=1./sqrt((cos(phi)/wx).^2+(sin(phi)/(wy)).^2);
        % r0=w0_phi/sqrt(2).*sqrt(-log(1 - rand_init(1,:)));
        px=wx/2 * randn(1,Np);
        py=wy/2 * randn(1,Np);
        r0=sqrt(px.^2+py.^2);
        fl = -w0/div_rad;
        PhStr.theta = -r0./fl;
        PhStr.cosphi = px./r0;
        PhStr.sinphi = py./r0;
        
%         PhStr.PreLoc(:,1) = px';
%         PhStr.PreLoc(:,2) = py';
%         PhStr.PreLoc(:,3) = pz;
        % Calculate next photon position x', y', z'  并转化坐标系
ux = sin(PhStr.theta).*PhStr.cosphi; %行向量
uy = sin(PhStr.theta).*PhStr.sinphi;
uz = cos(PhStr.theta);
%         uxyz_1 = [ux',uy',uz']*T01;
%         norm_uxyz_1 = sqrt(uxyz_1(:,1).^2 + uxyz_1(:,2).^2 +uxyz_1(:,3).^2 );
%         ux = (uxyz_1(:,1))';        uy=(uxyz_1(:,2))';      uz = ( uxyz_1(:,3) )'; 
        
%         scattering_events = ceil((c-a)*L*10);

%%发射机坐标系到RIS平面坐标系转换
Tx=0;Ty=La*sin(theta0);Tz=La*cos(theta0);
RIS_px=px+Tx;
RIS_py=py*cos(theta0-pi)+Ty;
RIS_pz=-py*sin(theta0-pi)+Tz;
RIS_ux=ux;
RIS_uy=uy*cos(theta0-pi)+uz*sin(theta0-pi);
RIS_uz=-uy*sin(theta0-pi)+uz*cos(theta0-pi);
ux=RIS_ux;
uy=RIS_uy;
uz=RIS_uz;
        PhStr.PreLoc(:,1) = RIS_px';
        PhStr.PreLoc(:,2) = RIS_py';
        PhStr.PreLoc(:,3) = RIS_pz';
        %% Interaction of photons with particles, see Figure 4.1
        while (1)
            % Photons still propogating?
            rand_seed;
            rand_array = rand(3,length(curr_set));
            % Choose random path length for traveling packets
            chi_dt = rand_array(1,:);
            r = -log(chi_dt)./c;
            
            Dist(curr_set) = Dist(curr_set) + r;
    
            % Choose scattring polar angles theta from CDF and phi
%             R1 = rand_array(2,:);%提取rand_array第二行所有列
%             R2 = rand_array(3,:);%提取rand_array第三行所有列
%             for i1 = 1:length(curr_set)
%                 ind = find(R1(i1)>=cdf_out);     %  已知 (x1,y1), (x2,y2), 以及yt, 求 xt， 进行插值处理
%                 if ind(end)+1<=length(cdf_out)
%                     x1 = rad(ind(end));
%                     y1 = cdf_out(ind(end)) ;
%                     x2 = rad(ind(end)+1);
%                     y2 = cdf_out(ind(end)+1);
%                     yt = R1(i1);
%                     xt = (yt-y1)*(x2-x1)/(y2-y1) + x1;
%                 else
%                     xt = rad(length(cdf_out));
%                 end
%                 PhStr.theta(curr_set(i1)) = xt;
%                 PhStr.phi(curr_set(i1)) = 2*pi*R2(i1);
%             end
    PhStr.theta(curr_set) = Hen_Green(length(curr_set),2); %1代表大气传输，2代表水中传输
    PhStr.phi(curr_set) = 2*pi*rand(1,length(curr_set));  
 
                
            if ~isempty(find(PhStr.theta(curr_set)>pi/2))
                tmp=1;
            end
            %              test_HG_CDF(curr_set,R1,R2,PhStr,cdf_out,rad);   %        测试HG散射相函数直接求theta的方法，并与累积分布函数方法比较
            
            x_up(curr_set) = r.*ux(curr_set);  % x'
        y_up(curr_set) = r.*uy(curr_set);  % y'
        z_up(curr_set) = r.*uz(curr_set);  % z'        
        
        % Set new position to x', y', z'
        PhStr.PostLoc(curr_set,:) = PhStr.PreLoc(curr_set,:) + [x_up(curr_set)', y_up(curr_set)', z_up(curr_set)'];
%         PhStr.PostLoc(curr_set,1) = PhStr.PreLoc(curr_set,1) + x_up(curr_set)';
%         PhStr.PostLoc(curr_set,2) = PhStr.PreLoc(curr_set,2) + y_up(curr_set)';
%         PhStr.PostLoc(curr_set,3) = PhStr.PreLoc(curr_set,3) + z_up(curr_set)';
        
        % Draw 2D-figure
%         cr = sqrt(PhStr.PostLoc(:,1).^2 + PhStr.PostLoc(:,2).^2 );
%         plot(cr ,PhStr.PostLoc(:,3),'k*'); hold on
%         vp = 50;  % view point
%         axis([-vp, vp,  -2*L, 2*L])
%         vec =[-vp:vp];
%         plot(vec,L*ones(size(vec)),'r-');
%         plot(vec,0*ones(size(vec)),'r-');
%         hold off
        
        % Photon in front of the Rx plane?
        % Choose photon packets which are still traveling, update curr_set
%         curr_set = find(PhStr.PostLoc(:,3)<L & PhStr.PostLoc(:,3)>=0) ;   % 查看累积的 z 方向距离是否超过 L, 且大于发送时的位置

        curr_set = find(PhStr.PostLoc(:,3)>0 & PhStr.PostLoc(:,3)<=Tz & PhStr.Wpost'~=0) ;   % 查看累积的 z 方向距离是否超过 L
%         curr_set = find(PhStr.PostLoc(:,3)<L  & PhStr.Wpost'~=0) ;   % 查看累积的 z 方向距离是否超过 L
        
        % Increment scatter number for survival photons
        PhStr.ScatNum(curr_set) = PhStr.ScatNum(curr_set) + 1;
        
        % Update PhStr.PreLoc(curr_set,:), 也即 curr_set 的PrePos位置上的光子仍然在传播，而其它位置上不需要更新，保持最后一次穿过Rx平面时候的位置
        PhStr.PreLoc(curr_set,:) = PhStr.PostLoc(curr_set,:);
        
        % Updata weight
        PhStr.Wpost(curr_set) = PhStr.Wpre(curr_set)*(1-a/c);
        
        % Rouletting

        indx_Th = find(PhStr.Wpost(curr_set)<Th );
        indx_other = setdiff(curr_set,curr_set(indx_Th));
        Roul_ind = [];
        if ~isempty(indx_Th)
            ind = curr_set(indx_Th);
            PhStr.Wpost(ind) = (rand(1,length(ind))<1/ah).*ah.*PhStr.Wpost(ind);
            Roul_ind = find(PhStr.Wpost(ind)~=0);
        end
        curr_set = [indx_other, curr_set(indx_Th(Roul_ind))];
        PhStr.Wpre = PhStr.Wpost;        
        
        % Photon packets in the Rx plane
        Rx_set =  setdiff([1:Np],curr_set);
        
        if isempty(curr_set)
            break;
        end
        

        % Set pre direction cosines
        rest_id = curr_set;
        
        id = find(abs(uz(curr_set))==1);
        if ~isempty(id)
            curr_id = curr_set(id);  % 注意，此处 id 对应的编号是 curr_set 中的序号，是个相对序号，所以要还原成绝对序号  curr_id
%             ux(curr_id) = sign(uz(curr_id)).*sin(PhStr.theta(curr_id)).*cos(PhStr.phi(curr_id));
%             uy(curr_id) = sign(uz(curr_id)).*sin(PhStr.theta(curr_id)).*sin(PhStr.phi(curr_id));
%             uz(curr_id) = sign(uz(curr_id)).*cos(PhStr.theta(curr_id));


            ux(curr_id) = sin(PhStr.theta(curr_id)).*cos(PhStr.phi(curr_id));
            uy(curr_id) = -sin(PhStr.theta(curr_id)).*sin(PhStr.phi(curr_id));
            uz(curr_id) = cos(PhStr.theta(curr_id));
            rest_id = setdiff(curr_set,curr_id);
        end
        % Calculate next photon position x', y', z', original
%         ux_pre = ux(rest_id);  uy_pre = uy(rest_id);  uz_pre = uz(rest_id);
%         us = cos(PhStr.theta(rest_id));
%         ux(rest_id) =sqrt(1-us.^2)./sqrt(1-uz_pre.^2).*   ( ux_pre.*uz_pre.*cos(PhStr.phi(rest_id)) - uy_pre.*sin(PhStr.phi(rest_id))  ) +...   %  ux(m+1) = sin(theta)*( ux(m)*uz(m)*cos(theta) - uy(m)*sin(phi) )/(sqrt(1-uz(m)^2)) + ux(m)*cos(theta)
%             ux_pre .*us;  
%         uy(rest_id) = sqrt(1-us.^2)./sqrt(1-uz_pre.^2).*  ( uy_pre.*uz_pre.*cos(PhStr.phi(rest_id)) +ux_pre.*sin(PhStr.phi(rest_id))  ) +...    %  uy(m+1) = sin(theta)*( uy(m)*uz(m)*cos(theta) + ux(m)*sin(phi) )/(sqrt(1-uz(m)^2)) + uy(m)*cos(theta)
%             uy_pre .*us;
%         uz(rest_id) =  -sqrt(1-us.^2).* cos(PhStr.phi(rest_id)) .* sqrt(1-uz_pre.^2) +...                                                       %  uz(m+1) = -sin(theta)*cos(phi)*(sqrt(1-uz(m)^2)) + uz(m)*cos(theta)
%             uz_pre .*us;
        
        
        %% new postion calculations
        ux_pre = ux(rest_id);  uy_pre = uy(rest_id);  uz_pre = uz(rest_id);
        us = cos(PhStr.theta(rest_id));
        ux(rest_id) =sqrt(1-us.^2)./sqrt(1-uz_pre.^2).*   ( uy_pre.*cos(PhStr.phi(rest_id)) + ux_pre.* uz_pre.*sin(PhStr.phi(rest_id))  ) +...   %  ux(m+1) = sin(theta)*( uy(m)*cos(theta) + ux(m)*uz(m)*sin(phi) )/(sqrt(1-uz(m)^2)) + ux(m)*cos(theta)
            ux_pre .*us;
        uy(rest_id) = sqrt(1-us.^2)./sqrt(1-uz_pre.^2).*  ( -ux_pre.*cos(PhStr.phi(rest_id)) +uy_pre.* uz_pre.*sin(PhStr.phi(rest_id))  ) +...   %  uy(m+1) = sin(theta)*( -ux(m)*cos(theta) + uy(m)*uz(m)*sin(phi) )/(sqrt(1-uz(m)^2)) + ux(m)*cos(theta)
            uy_pre .*us;
        uz(rest_id) =  -sqrt(1-us.^2).* sin(PhStr.phi(rest_id)) .* sqrt(1-uz_pre.^2) +...                                                        %  uz(m+1) = -sin(theta)*sin(phi)*(sqrt(1-uz(m)^2)) + uz(m)*cos(theta)
            uz_pre .*us;

        if ~isempty(find((-uy*sin(pi-theta0)+uz*cos(pi-theta0))<0))
            tmp=1;
        end

        
        
        count = count +1;
        
        % Draw 3D-figure
        %             plot3(PhStr.PostLoc(:,1),PhStr.PostLoc(:,2),PhStr.PostLoc(:,3),'k*'); hold on
        %             vp = 100;  % view point
        %             axis([-vp, vp,  -vp, vp, 0, 2*L])
        %             draw_plane(L);  hold off
        
    end
        
        %% Goto Receiver
        % 方法一：流程为先求空间直线与Rx平面的交点；然后根据交点位置确定离Rx接收平面原点的距离
        % 空间直线与Rx平面的交点，http://blog.csdn.net/abcjennifer/article/details/6688080
        % 设穿过Rx平面后的点坐标为(x_up,y_up,z_up),对应PhStr.PostLoc； 散射前的坐标为(x,y,z),对应PhStr.PreLoc 求这两点确定的空间直线与平面的交点(xp,yp, zp)
        % 1) 空间直线的参数方程可表示为：
        % xp = x + v1*t;
        % yp = y + v2*t;
        % zp = z + v3*t
        % 其中 v1 = (x_up - x),  v2 = (y_up - y), v3 = (z_up - z)
        % 2) 平面点法式方程，平面上的任意一点( x_up,y_up,L), 法向量 (0, 0, z_up-L)
        % 0*(xp - x_up ) + 0*(yp- y_up) + (z_up-L)*(zp - L) = 0
        % 3) 实际上得出 zp = L, 从而 t =  (L-z)/v3 = (L - z)/(z_up - z),  并进一步求得 xp, yp
        
        
        
    %     %     uo = find(PhStr.PostLoc(:,3)<0);
    %     ta = ( L-PhStr.PreLoc(:,3) )./(PhStr.PostLoc(:,3) - PhStr.PreLoc(:,3));
    %     xa = PhStr.PreLoc(:,1) + ( PhStr.PostLoc(:,1) - PhStr.PreLoc(:,1)).*ta;
    %     ya = PhStr.PreLoc(:,2) + ( PhStr.PostLoc(:,2) - PhStr.PreLoc(:,2)).*ta;
    %     %     RALL = sqrt(xa.^2 + ya.^2);
    %     
    %     Updis = PhStr.PostLoc(:,3) - L;
    %     Eucdis = sqrt(  (PhStr.PostLoc(:,1) -  xa).^2 + (PhStr.PostLoc(:,2) - ya).^2  + (PhStr.PostLoc(:,3) - L).^2  ) ;
    %     incid_rad = acos(Updis./Eucdis);
    %     incid_ang = incid_rad/pi*180;
    %     
    %     ui = find(PhStr.PostLoc(:,3)>0 & incid_ang<=FOV);  % 注意：只统计到达 zp=L 平面的所有光子，而小于zp=0 的不统计, 且满足FOV要求
    %     
    %     t = ( L-PhStr.PreLoc(ui,3) )./(PhStr.PostLoc(ui,3) - PhStr.PreLoc(ui,3));
    %     xp = PhStr.PreLoc(ui,1) + ( PhStr.PostLoc(ui,1) - PhStr.PreLoc(ui,1)).*t;
    %     yp = PhStr.PreLoc(ui,2) + ( PhStr.PostLoc(ui,2) - PhStr.PreLoc(ui,2)).*t;
    %     zp = L;
    %     Rd = sqrt(xp.^2 + yp.^2);
     
    
        % 方法二：
        uxyz_3=[ux',uy',uz'];   %海面坐标系下方向向量
        ux_3 = uxyz_3(:,1);     uy_3 = uxyz_3(:,2);     uz_3 = uxyz_3(:,3);
        ui = find(PhStr.PostLoc(:,3)<0 );  % 注意：只统计到达 zp=L 平面的所有光子，而小于zp=0 的不统计, 且满足FOV要求    
        zd = 0 - PhStr.PreLoc(ui,3);
        yd = zd.*(uy(ui)./uz(ui))';
        xd = zd.*(ux(ui)./uz(ui))';
        rec_x = xd +  PhStr.PreLoc(ui,1);
        rec_y = yd +  PhStr.PreLoc(ui,2);
        rec_w = PhStr.Wpost(ui);
        rec_ux = ux_3(ui);
        rec_uy = uy_3(ui);
        rec_uz = uz_3(ui);
    
    
%         xyz_3 = [PhStr.PostLoc(:,1),PhStr.PostLoc(:,2),PhStr.PostLoc(:,3)]; %海面坐标系下位置坐标
%         uxyz_3=[RIS_ux',RIS_uy',RIS_uz'];   %海面坐标系下方向向量
%         ux_3 = uxyz_3(:,1);     uy_3 = uxyz_3(:,2);     uz_3 = uxyz_3(:,3);
%         ui = find(xyz_3(:,3)<0);   
%         zd  = 0 - xyz_3(ui,3);%大于0
%         yd  = zd.*(uy_3(ui)./uz_3(ui));%zd/uz 小于0
%         xd  = zd.*(ux_3(ui)./uz_3(ui));
%         rec_x  = xd + xyz_3(ui,1);  %海平面坐标系下坐标
%         rec_y  = yd + xyz_3(ui,2); %
%         rec_w = PhStr.Wpost(ui);
%         rec_ux = ux_3(ui);
%         rec_uy = uy_3(ui);
%         rec_uz = uz_3(ui);
        
%         L=0;
%         ui = find(PhStr.PostLoc(:,3)>L);  % 注意：只统计到达 zp=L 平面的所有光子，而小于zp=0 的不统计, 且满足FOV要求     绝对序号
%         ui = find(PhStr.PostLoc(:,3)<L);
        %% test onlt
%         z_all_post = PhStr.PostLoc(:,3)';
%         z_all_pre  = PhStr.PreLoc(:,3)';
%         z_initial = PhStr.PostLoc(ui,3)';
%         z_initial_pre = PhStr.PreLoc(ui,3)';
        %%  这部分为到达目标平面的光子
%         zd  = L - PhStr.PostLoc(ui,3); %列向量   此场景下大于0
%         yd  = zd.*(uy(ui)./uz(ui))';
%         xd  = zd.*(ux(ui)./uz(ui))';
%         rec_x  = xd + PhStr.PostLoc(ui,1);
%         rec_y  = yd + PhStr.PostLoc(ui,2); %
%         rec_w = PhStr.Wpost(ui);
%         rec_ux = ux(ui);
%         rec_uy = uy(ui);
%         rec_uz = uz(ui);
        para_1.x = rec_x;       para_1.y = rec_y;       para_1.z = zeros(size(rec_x));
        para_1.ux = rec_ux;     para_1.uy = rec_uy;     para_1.uz = rec_uz;
        para_1.w   = rec_w;
%        rec_rx = rec_x; 
%        rec_ry = rec_y; 
%        rec_rw = rec_w;
        %%
 rec_r = sqrt(rec_x.^2 + rec_y.^2); %% 到达接受面的光子的r            
 R0 = linspace(0,1,100);


for i = 1:length(R0)-1
s=pi*(R0(i+1)^2  - R0(i)^2 );
rx_id = find(rec_r>R0(i) & rec_r<=R0(i+1));
rx_pow(nn,i)=sum(PhStr.Wpost(ui(rx_id)))/(Np)/s;
end
%         Rd = rec_x.^2 + rec_y.^2;
%         %进入棱镜范围内的光子 ，输出其 坐标，方向向量、权重
%         T = find(Rd(:,1)<= R^2);%再次判断条件 到达接受面的光子中，哪些进入了半径范围内
%         rec_rx  = rec_x (T);
%         rec_ry  = rec_y (T);
%         rec_rw = rec_w(T);
%         r_ux    =  ux(ui(T));
%         r_uy    =  uy(ui(T));
%         r_uz    =  uz(ui(T)); 
%         
%         time(T) = Dist(T)/speed*1e9;
%         time_reach = sort(time(T));
%         Timereach = [Timereach,time_reach];
%         [~, sorted_indices] = sort(time_reach);
%         T_sorted = T(sorted_indices);
%         Wpost_reach = PhStr.Wpost(T_sorted);
%         Wreach = [Wreach,Wpost_reach];
%         arrived_index = ui(T); %绝对序号
% %         Dist(arrived_index) =  
%         arrived_time  = Dist(arrived_index)./speed*1e9;
%         arrived_w     = PhStr.Wpost(arrived_index);
%         [time_reach,index] = sort(arrived_time);
%         Timereach = [Timereach,time_reach];
%         Wpost_reach = arrived_w(index);
%         Wreach    = [Wreach,Wpost_reach];
% 
%         min_time_number = min(time_reach);
%         max_time_number = max(time_reach);
%         Time_number = min_time_number:0.1:max_time_number;
%         for it = 1:length(Time_number)
%             current_time_number = Time_number(it); % 当前时间
%             valid_indices_number = find(abs(time_reach - current_time_number) <= 0.02);
%             Wrow{nn}(it) = sum(Wpost_reach(valid_indices_number))/(Np*P0);
%         end
%     end
% 
%         [Timereach , initial_index] = sort(Timereach);
%         Wreach = Wreach(initial_index);
    end
rx_pow_out=mean(rx_pow,1);
tmp=1;