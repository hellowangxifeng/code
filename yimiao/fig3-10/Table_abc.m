function [a,b,c] = Table_abc(water_type)
%% water type, see TABLE I
c3_2    =	[2.6, 3.2];
c2_675 =   [2.1734 , 2.675];
c1_92   =	[1.56, 1.92];
c1_475 =    [1.1984, 1.475];
c1_0 =      [0.8125, 1];
c0_5 =      [0.406, 0.5];
c0_275   =	[0.2234, 0.275];
c0_15  = [0.08,0.15];%远洋海水的散射和衰减系数
c0_305  = [0.216,0.305];%近岸海水的散射和衰减系数
c0_151 = [0.037,0.151];
c0_398 = [0.219,0.398];
c2_190  = [1.824, 2.190];
c2_17 = [1.875,2.17];
if strcmp(water_type,'c3.2')
    paras = c3_2;  % 选用 phase scattering function of hankel transform:  ph_Pet, ph_HG, ph_FF
elseif strcmp(water_type,'c2.675')
    paras = c2_675;
elseif strcmp(water_type,'c1.92')
    paras = c1_92;
elseif strcmp(water_type,'c1.475')
    paras = c1_475;
elseif strcmp(water_type,'c1.0')
    paras = c1_0;
elseif strcmp(water_type,'c0.5')
    paras = c0_5;
elseif strcmp(water_type,'c0.275')
    paras = c0_275;
elseif strcmp(water_type,'c0.15')
    paras = c0_15;
elseif strcmp(water_type,'c0.305')
    paras = c0_305;
elseif strcmp(water_type,'c2.17')
    paras = c2_17;
elseif strcmp(water_type,'c0.151')
    paras = c0_151;
elseif strcmp(water_type,'c0.398')
    paras = c0_398;
elseif strcmp(water_type,'c2.190')
    paras = c2_190;

end

b = paras(1);
c = paras(2);
a = c - b;