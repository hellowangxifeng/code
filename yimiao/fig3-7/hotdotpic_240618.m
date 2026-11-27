function w= hotdotpic_240618(x,y,rx,ry,rw,Np,dmax,m1_g)
%HOTDOTPIC 此处显示有关此函数的摘要
%   此处显示详细说明
%m
% len = 100;
% h = hotdotpic(x_dir,rec_x',rec_y', rec_w,1e7,100); 
% y=x; 
len = length(x) -1 ;
w=zeros(len,len);
%(i,j)
for i=1: len
    for j = 1: len
           i,j
           index = find(   rx>x(i)  &  rx<=x(i+1) &  ry>y(j)  &  ry<=y(j+1) );
%            w(i,j)  = sum(rw(index))/Np;
           w(i,j)  = sum(rw(index))/Np/((2*dmax/m1_g)^2);
    end
end
% h=heatmap(w);

end

