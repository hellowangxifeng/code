% 2013_Monte-Carlo-Based Channel Characterization for Underwater Optical Communication Systems
% function [MC_BSF, Vz, RX_POW] = Photon_Propagation_Parallel(P0,number, water_type,psf_type,L, Rx_r, dr, div_rad, lambda)
function [aaa] = MC_sea_in_3axis_250319_2(P0, number, water_type,L ,x,y,w,ux,uy,uz,fi0)


% mode =0 用上一阶段大气结束时的光子方向向量 
% mode=1 根据光子的位置x和y重新生成方向向量



% psf_type='HG';
%% generate CDF of phase scattering function\\分布函数
index=2; %1代表大气g=0.86     2代表水下 g=0.924
% [cdf_out, rad]= psf_cdf(psf_type,index);
% loglog(rad,cdf_out,'k-','LineWidth',2)   % 对比 simulation, modeling, and design of underwater optical communication systems, Fig. 4.27

[a,b,c] = Table_abc(water_type);
% max_time = 0;
% min_time = 0;
speed = 2.25*10^8;  % speed of light in water %m/s
% R=Rw;
% Timereach=[];
% Wreach=[];
% Wsum=[];
% Wrow=cell(number,1);
% Wsample=zeros(number,sample_point);
% rx_pow = zeros(number,length(Rx_r));
% mc_bsf = zeros(number, length(Rx_r));
trials = 1;
for test=1:trials
    for nn=1:number
        %% Initialization
%         nn
%         Np = 1;
        Np = length(w);  % number of photons
        FOV = 180;
    %     div_deg = 1;  %0.01;  % divergence angle in degree
    %    div_rad = div_deg/180*pi; % divergence angle in radians
    %     lambda = 532e-9; % wavelength
    
        PreLoc = zeros(Np,3);  % pre-position
        PreLoc(: ,1)  = x; PreLoc(:,2) = y ;
        PostLoc = zeros(Np,3);   % initial positions for (x,y,z)
        PostLoc(:,1) = x; PostLoc(:,2) = y ;
        Wpre = w;  % pre weight, Initially, each photon is launched into the medium with unit weight. 
        Wpost = w;  % post weight,  Initially, each photon is launched into the medium with unit weight.
        time = zeros(1,Np);
        Dist = zeros(1,Np);    % Photon's distance traveled
        theta =zeros(1,Np); % polar angle in radians
        phi = zeros(1,Np);  % azimuthal angle in radians
        ScatNum  = zeros(1,Np);  % count scatter number for each photon
        
        PhStr = struct('Wpre',Wpre,'Wpost',Wpost, 'PreLoc', PreLoc, 'PostLoc',PostLoc,...
            'Dist', Dist, 'theta',theta,'phi',phi, 'ScatNum', ScatNum);  % photon characteristics
        
%         Th = 1e-6;   % weight dropping threshold
        Th = 1e-4;
        ah = 10; % rouletting threshold
        count = 1;
        curr_set = [1:Np];
            
        % Initiate scattering angles
%         w0 =  lambda/(pi*div_rad);  % refer to :  https://en.wikipedia.org/wiki/Gaussian_beam#Radius_of_curvature,    Beam divergence
%     
        rand_seed;
        rand_init = rand(2,Np);
%         wz = w0*(1+(lambda*L/(pi*w0^2))^2)^(1/2);
%         Vz = wz^2/4;
%         r0 = w0/sqrt(2).*sqrt(-log(1 - rand_init(1,:))); % array of selected values r1
%         fl = -w0/div_rad;
%   
%         PhStr.theta = -r0./fl;    % polar angle in radians
%         PhStr.phi = 2*pi*rand_init(2,:);        % azimuthal angle in radians
         
%           fl = Wz/div_rad;

%         px=r0'.*cos(PhStr.phi)';
%         py=r0'.*sin(PhStr.phi)';
%         PhStr.PreLoc(:,1) = px;
%         PhStr.PreLoc(:,2) = py;
        
        % Calculate next photon position x', y', z'
%         ux,uy,uz 由api接口给出
% mode =0 用上一阶段大气结束时的光子方向向量 
% mode=1 根据光子的位置x和y重新生成方向向量
%     if mode == 1
%             r0=sqrt(x.^2 + y.^2);
%            PhStr.theta =  (r0/fl)';
%            PhStr.phi = 2*pi*rand_init(2,:);  
%         ux = sin(PhStr.theta).*cos(PhStr.phi);
%         uy = sin(PhStr.theta).*sin(PhStr.phi);
%         uz = cos(PhStr.theta);
%     end
        scattering_events = ceil((c-a)*L*10);
        %% Interaction of photons with particles, see Figure 4.1
        while (1)
            % Photons still propogating?
            rand_seed;
            rand_array = rand(3,length(curr_set));
            % Choose random path length for traveling packets
            chi_dt = rand_array(1,:);
            r = -log(chi_dt)./c;
            
            Dist(curr_set) = Dist(curr_set) + r;
    
%             Choose scattring polar angles theta from CDF and phi
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
%     PhStr.theta(curr_set) = Hen_Green(length(curr_set),2);
%     PhStr.phi(curr_set) = 2*pi*rand(1,length(curr_set));  
    PhStr.theta(curr_set) = Hen_Green(length(curr_set),2); %1代表大气传输，2代表水中传输
    PhStr.phi(curr_set) = 2*pi*rand(1,length(curr_set));
                
            if ~isempty(find(PhStr.theta(curr_set)>pi/2))
                tmp=1;
            end
            %              test_HG_CDF(curr_set,R1,R2,PhStr,cdf_out,rad);   %        测试HG散射相函数直接求theta的方法，并与累积分布函数方法比较
            
            x_up(curr_set) = r'.*ux(curr_set)';  % x'
            y_up(curr_set) = r'.*uy(curr_set)';  % y'
            z_up(curr_set) = r'.*uz(curr_set)';  % z'        
            
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
    
            curr_set = find(PhStr.PostLoc(:,3)<L & (-sin(fi0)*PhStr.PostLoc(:,1)+cos(fi0)*PhStr.PostLoc(:,3))>=0 & PhStr.Wpost'~=0) ;   % 查看累积的 z 方向距离是否超过 L
    %         curr_set = find(PhStr.PostLoc(:,3)<L  & PhStr.Wpost'~=0) ;   % 查看累积的 z 方向距离是否超过 L
            
            % Increment scatter number for survival photons
            PhStr.ScatNum(curr_set) = PhStr.ScatNum(curr_set) + 1;
            
            % Update PhStr.PreLoc(curr_set,:), 也即 curr_set 的PrePos位置上的光子仍然在传播，而其它位置上不需要更新，保持最后一次穿过Rx平面时候的位置
            PhStr.PreLoc(curr_set,:) = PhStr.PostLoc(curr_set,:);
            
            % Update weight
            PhStr.Wpost(curr_set) = PhStr.Wpre(curr_set)*(1-a/c);
            
            % Rouletting
    
            indx_Th = find(PhStr.Wpost(curr_set)<Th );%权值小于阈值的序列号
            indx_other = setdiff(curr_set,curr_set(indx_Th));%找出curr_set中权值大于阈值的序列号
            Roul_ind = [];
            if ~isempty(indx_Th)

%                         length(curr_set);

                ind = curr_set(indx_Th);
                PhStr.Wpost(ind) = (rand(1,length(ind))<1/ah).*ah.*PhStr.Wpost(ind);
                Roul_ind = find(PhStr.Wpost(ind)~=0);
            end
            curr_set = [indx_other; curr_set(indx_Th(Roul_ind))];
            PhStr.Wpre = PhStr.Wpost; %权值更新       
            
            A=PhStr.Wpost;
            Z=PhStr.PostLoc(:,3)';
            % Photon packets in the Rx plane
            Rx_set =  setdiff([1:Np],curr_set);%平面接受光子序号
            
           
%                  length(curr_set)
        
            
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
        
            if ~isempty(find(uz<0))
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
            ui = find(PhStr.PostLoc(:,3)>L);  % 注意：只统计到达 zp=L 平面的所有光子，而小于zp=0 的不统计, 且满足FOV要求     绝对序号
            %修正
            zd = L - PhStr.PostLoc(ui,3); %列向量
            yd = zd.*(uy(ui)./uz(ui))';
            xd = zd.*(ux(ui)./uz(ui))';
            rec_x = xd +  PhStr.PostLoc(ui,1);
            rec_y = yd +  PhStr.PostLoc(ui,2); %
            
            PhStr.PostLoc(ui,1) = rec_x;
            PhStr.PostLoc(ui,2) = rec_y;
            Dist(ui) =  Dist(ui) + (zd./uz(ui)')';
            
            %输出
            rec_w  = PhStr.Wpost(ui);
           
% %             time   = rt(ui) +Dist(ui)./speed*1e9;
%             
%             Rd = (rec_x - 0).^2 + (rec_y -0) .^2;
%             Rd = Rd';
%             T = find(Rd<= R^2);%再次判断条件
%             
%             aaa.x = rec_x(T);
%             aaa.y = rec_y(T);
%             aaa.w = rec_w(T);
%             aaa.t = time(T);
%             aaa.ux = ux(ui(T));
%             aaa.uy = uy(ui(T));
%             aaa.uz = uz(ui(T));
%              
%                rec_r = sqrt(rec_x.^2 + rec_y.^2); %% 到达接受面的光子的r
                
%                R0 = linspace(0,1,500);
%             rx_id=find(rec_x <= R0(i)+r_a & rec_x >= R0(i)-r_a & rec_y.^2+(rec_x-R0(i)).^2<=r_a.^2);
%             for i = 1:length(R0)-1
%                 s=pi*(R0(i+1)^2  - R0(i)^2 );
%                 rx_id = find(rec_r>R0(i) & rec_r<=R0(i+1));
%                 rx_pow(nn,i)=sum(PhStr.Wpost(ui(rx_id)))/(1e7)/s;
%             end         
            aaa.x=rec_x;
            aaa.y=rec_y;
            aaa.w=rec_w;

    end
end
% RX_POW = mean(rx_pow,1);
% % Ir = 10*log10(Rx_pow/Np);
% MC_BSF = mean(mc_bsf,1);
% bsf_r = mean(rx_pow,1);
tmp=1;