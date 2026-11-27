function ellipse=rp_spot_ellipse(spot)
%RP_SPOT_ELLIPSE Second-moment waist ellipse, including its actual orientation.
arguments
    spot (1,1) struct
end
C=(spot.covarianceM2+spot.covarianceM2.')/2;
[vectors,values]=eig(C,'vector');[values,order]=sort(values,'descend');
assert(all(values>=-1e-12),'rp:SpotCovariance','Spot covariance is not positive semidefinite.');
radii=2*sqrt(max(values,0)); vectors=vectors(:,order);
angle=linspace(0,2*pi,361);
points=(vectors*(radii.*[cos(angle);sin(angle)])).';
ellipse=struct('centerXY',spot.centerXY,'semiAxesM',radii.', ...
    'directions',vectors,'offsetXY',points,'axisEndsXY',vectors*diag(radii));
end
