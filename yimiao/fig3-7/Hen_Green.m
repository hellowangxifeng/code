function angle = Hen_Green(length,judge)
if judge == 1
    g = 0.86;
else
    g = 0.924;
end
val= (1+g^2-((1-g^2)./(1-g+2*g*rand(1,length))).^2)*(1/(2*g));
angle=acos(val);
end


