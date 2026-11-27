function bound=rp_nonstationary_rect_bound(amplitude,gx,gy,x,y)
%RP_NONSTATIONARY_RECT_BOUND Two integrations by parts, rectangular aperture.
% For integral A exp(i Phi), gx=Phi_x and gy=Phi_y must not vanish.
assert(all(abs(gx)>0 & abs(gy)>0,'all'),'rp:StationaryPoint','Phase derivative vanishes.');
dx=x(2)-x(1); dy=y(2)-y(1);
corner=amplitude./gx./gy;
[~,edgeDerivative]=gradient(corner,dx,dy);
[firstDerivative,~]=gradient(amplitude./gx,dx,dy);
remainder=firstDerivative./gy;
[~,mixedDerivative]=gradient(remainder,dx,dy);
bound=sum(abs(corner([1 end],[1 end])),'all')+ ...
    trapz(y,abs(edgeDerivative(:,1))+abs(edgeDerivative(:,end)))+ ...
    trapz(x,abs(remainder(1,:))+abs(remainder(end,:)))+ ...
    trapz(y,trapz(x,abs(mixedDerivative),2));
end
