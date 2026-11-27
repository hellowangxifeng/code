function [rx_pow] = getbsf_in_fi_240617(aaa,fi,n)
%GETBSF_IN_FI 此处显示有关此函数的摘要
r=linspace(0,0.5,101);  %2cm×2cm的网格为极限0.01
cell_a=0.02;
x=r*cos(fi);y=r*sin(fi);
for i=1:length(r)
%     i
        index  =  find(  abs( aaa.x -x(i) )<=cell_a/2 & abs( aaa.y -y(i) )<=cell_a/2  );
        rx_pow(i) = sum(aaa.w(index))/n/(cell_a^2);
end

end

