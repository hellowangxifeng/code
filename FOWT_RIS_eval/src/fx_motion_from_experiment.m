function m = fx_motion_from_experiment(file,opts)
%FX_MOTION_FROM_EXPERIMENT Build an RIS_float-compatible motion struct from a
%   measured Qualisys 6-DOF record.
%
%   m = FX_MOTION_FROM_EXPERIMENT(file) reads a FOWTmotion_*.mat file and
%   returns the minimal rigid-body struct used by rp_geometry/rp_snapshot:
%       t, rotationMatrix, bottomCenter, normal, boardCenter, centerOfMass,
%       rollDeg/pitchDeg/yawDeg, heave/surge/sway, parameters
%
%   The record is centred on its mean horizontal position and mean heave, so
%   the mean pose is the design pose. This is the one-time trim calibration
%   discussed in the paper. Pass opts.center=true to keep absolute positions.
%
%   opts.resampleHz (default 100), opts.thickness (default 0.01),
%   opts.z0 (default -0.0334, the RIS mean bottom height), opts.name.
arguments
    file (1,:) char {mustBeFile}
    opts.resampleHz (1,1) double {mustBePositive} = 100
    opts.thickness (1,1) double {mustBePositive} = 0.01
    opts.z0 (1,1) double = -0.0334
    opts.center (1,1) logical = true
    opts.name (1,:) char = ''
end
S=load(file);
names=fieldnames(S);
names=names(~startsWith(names,'__'));
assert(numel(names)==1,'fx:Motion','Expected one Measurement struct in %s.',file);
M=S.(names{1});
rb=M.RigidBodies;
pos=double(squeeze(rb.Positions)).';      % N x 3, mm
rpy=double(squeeze(rb.RPYs)).';           % N x 3, deg
fs=double(M.FrameRate);
n=min([size(pos,1),size(rpy,1),double(M.Frames)]);
pos=pos(1:n,:); rpy=rpy(1:n,:);
% linear NaN fill
for j=1:3
    pos(:,j)=fillOne(pos(:,j));
    rpy(:,j)=fillOne(rpy(:,j));
end
t=(0:n-1).'/fs;
% resample to a uniform grid
tq=(0:1/opts.resampleHz:t(end)).';
pos=interp1(t,pos,tq,'linear','extrap');
rpy=interp1(t,rpy,tq,'linear','extrap');
t=tq; n=numel(t);
if opts.center
    pos(:,1)=pos(:,1)-mean(pos(:,1));
    pos(:,2)=pos(:,2)-mean(pos(:,2));
    pos(:,3)=pos(:,3)-mean(pos(:,3));
end
RRis=zeros(3,3,n); normal=zeros(3,n); bottom=zeros(3,n); board=zeros(3,n); com=zeros(3,n);
for k=1:n
    R=rotationZYX(deg2rad(rpy(k,:).'));
    RRis(:,:,k)=R;
    normal(:,k)=R(:,3);
    bottom(:,k)=[pos(k,1)*1e-3;pos(k,2)*1e-3;opts.z0+pos(k,3)*1e-3];
    board(:,k)=bottom(:,k)+opts.thickness*normal(:,k);
    com(:,k)=bottom(:,k)+0.5*opts.thickness*normal(:,k);
end
p=get_default_RIS_params();
m=struct('t',t,'rotationMatrix',RRis,'bottomCenter',bottom, ...
    'normal',normal,'boardCenter',board,'centerOfMass',com, ...
    'rollDeg',rpy(:,1),'pitchDeg',rpy(:,2),'yawDeg',rpy(:,3), ...
    'heave',pos(:,3)*1e-3,'surge',pos(:,1)*1e-3,'sway',pos(:,2)*1e-3, ...
    'parameters',p,'source',file,'name',opts.name, ...
    'note','measured Qualisys 6-DOF, centred; not the 1 m RIS hydrodynamics');
end

function y=fillOne(y)
bad=~isfinite(y);
if ~any(bad), return; end
idx=(1:numel(y)).';
y(bad)=interp1(idx(~bad),y(~bad),idx(bad),'linear','extrap');
end

function R=rotationZYX(e)
% Body-to-world Rz(yaw) Ry(pitch) Rx(roll); e=[roll;pitch;yaw].
r=e(1); p=e(2); y=e(3);
Rx=[1 0 0;0 cos(r) -sin(r);0 sin(r) cos(r)];
Ry=[cos(p) 0 sin(p);0 1 0;-sin(p) 0 cos(p)];
Rz=[cos(y) -sin(y) 0;sin(y) cos(y) 0;0 0 1];
R=Rz*Ry*Rx;
end
