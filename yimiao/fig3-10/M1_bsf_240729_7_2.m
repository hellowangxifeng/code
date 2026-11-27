function [out2] =M1_bsf_240729_7_2(P0,La,Lw,theta,fi,div_rad,water_type,x,y)

% x=linspace(-0.2,0.2,100);
[x1,y1]=meshgrid(x,y);
% %%get initial bsf
% %  La=100;     theta=30/180*pi;           iii=3;
% beta1 = beta;
% % div_rad=1e-3;
% % A=beta1(1); sigma1_1=abs(beta1(4)); sigma2_1 = abs(beta1(5));
%  A=beta1(1); sigma1_1=abs(beta1(2)); sigma2_1 = abs(beta1(3));
% d_left = La*sin(div_rad)/sin(theta-div_rad);
% sin_fi = 0.75*sin(pi/2 - theta); %入射角为pi/2 - theta
% L0=d_left*sin_fi;
% %%计算水下发散角
% theta_in = [pi/2-(theta + div_rad), pi/2 - theta,  pi/2-(theta - div_rad)]; %入射角
% % theta_in=[pi/2 - theta,pi/2 - theta,pi/2 - theta];=
% sin_fi_sea = 0.75*sin(theta_in);%三个折射角正弦
% 
% fi_sea = asin(sin_fi_sea);
% div_rad_sea_y = atan(sin(div_rad)/cos(div_rad+theta));  
% div_rad_sea_y = div_rad/cos(theta);  
% div_rad_sea_x =div_rad;
%% begin caculate bsf
% if water_type=='c0.398'
%     kk=1;
% elseif water_type=='c0.151'
%     kk=0;
% end
[~,b,c] = Table_abc(water_type);
k1 = -0.07973;
ek2 = 1.05;
% L=10; A=1; sigma1=1;sigma2=1;
% L=Lw - L0;
% L=Lw;
% sigma1 = sigma1_1+(div_rad)*(L-L0)/2; % wz = 2 sigma;  %%应该成比例放大
% sigma1 = sigma1_1+0.0008+sigma2_1/sigma1_1*(div_rad_sea_y)*(L-L0); % wz = 2 sigma;  %%应该成比例放大
lda=523e-9;
w0=lda/(pi*div_rad);
WL_x =@(d) w0*sqrt(1+(lda*d/(pi*w0^2)).^2); % wz = 2 sigma;  %%应该成比例放大
WL_y =@(d) w0*sqrt(1+(lda*d/(pi*w0^2)).^2);
% y1=exp(b*ek2*(exp(k1.*x.*L)-1)./(k1.*x));
% 发射机光源 参数 A sigma1，sigma2
E1=@(wx,wy,d,d2) P0*exp(-((WL_x(d)+div_rad*d2)^2*wx.^2 +(WL_y(d)*cos(fi)/cos(theta)+div_rad*d2)^2*wy.^2)/8);
%O                                  VF
S = @(wx,wy,d,d2) exp(-c*(d+d2)).*exp( b*ek2*exp(k1.*sqrt(wx.^2 + (wy).^2).*(d2)).*(exp(k1.*sqrt(wx.^2 + (wy*cos(fi)/cos(theta)).^2).*(d))-exp(k1.*sqrt(wx.^2 + (wy*cos(fi)/cos(theta)).^2).*0))./(k1.*sqrt(wx.^2 + (wy.*cos(fi)/cos(theta)).^2))).*exp(  b*ek2*(   exp(k1.*sqrt(wx.^2 + (wy).^2).*(d2))-exp(k1.*sqrt(wx.^2 + (wy).^2).*0))./(k1.*sqrt(wx.^2 + (wy).^2)  )  );
F=@(wx,wy,d,d2) 1/4/pi/pi*E1(wx,wy,d,d2).*S(wx,wy,d,d2);
%sinb=(L*sin(pi/2-theta)/(sqrt((y/cos(fi))^2+L^2-2*(y/cos(fi))*L*cos(pi/2-theta))));
A=La*cos(fi)/cos(theta);
B=-tan(theta);
C=A+Lw;
D=tan(fi)+B;
%d1=A+B/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y));
%d2=Lw+tan(fi)/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y));
bsf_end = @(x,y)  arrayfun(@(x,y)   integral2(@(wx,wy) (F(wx,wy,(La-sin(theta)/cos(fi)/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y))),(Lw+tan(fi)./(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y)))).*exp(1i*x.*wx+1i*y.*wy)*((cos(fi)/cos(theta))/(cos(theta-fi)+sin(theta-fi)*sqrt(1-(La*sin(pi/2-theta)/(sqrt((1/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y))/cos(fi))^2+La^2-2*(1/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y))/cos(fi))*La*cos(pi/2-theta))))^2)/(La*sin(pi/2-theta)/(sqrt((1/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y))/cos(fi))^2+La^2-2*(1/(2*D)*(B*y-C+sqrt(C^2+B^2*y.^2+(4*A*D-2*B*C)*y))/cos(fi))*La*cos(pi/2-theta))))))),-inf,inf,-inf,inf,'AbsTol' ,...
10^(-5) , 'RelTol' ,10^(-5) ),x,y);
% bsf_end = @(x,y)  arrayfun(@(x,y)   integral2(@(wx,wy) (F(wx,wy,(L-y.*sin(theta))).*exp(1i*x.*wx+1i*y.*wy*cos(theta))*(cos(theta))),-inf,inf,-inf,inf,'AbsTol' ,...
% 10^(-5) , 'RelTol' ,10^(-5) ),x,y);
%sinb=(-y*cos(theta)/(L-y*sin(theta)))
% bsf_end = @(x,y)  arrayfun(@(x,y)   integral2(@(wx,wy) (F(wx,wy,(L-y.*sin(theta))).*exp(1i*x.*wx+1i*y.*wy*cos(theta))*(cos(theta)-sin(theta)*(-y*cos(theta)/(L-y*sin(theta)))/sqrt(1-(-y*cos(theta)/(L-y*sin(theta)))^2))),-inf,inf,-inf,inf,'AbsTol' ,...
% 10^(-5) , 'RelTol' ,10^(-5) ),x,y);
%
x_row = reshape(x1,1,length(x1)*length(x1)); %将X按列展开为一行 [123;456] ---> [132536]
y_row = reshape(y1,1,length(y1)*length(y1));  %这里用x1 y1 相当于重新划分区域 与之前的X_M Y_M 无关

out=bsf_end(x_row,y_row); %985.6秒 约16分钟
%%
out2 = reshape(abs(out),size(x1)); 


end

