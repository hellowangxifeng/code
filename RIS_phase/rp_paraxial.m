function o = rp_paraxial(c,g,reference,targets)
%RP_PARAXIAL Exact central layered rays and canonical second-order OPL matrices.
arguments
    c (1,1) struct
    g (1,1) struct
    reference (:,2) double {mustBeReal,mustBeFinite}
    targets (:,2) double {mustBeReal,mustBeFinite}
end
n=max(size(reference,1),size(targets,1));
if size(reference,1)==1, reference=repmat(reference,n,1); end
if size(targets,1)==1, targets=repmat(targets,n,1); end
assert(size(reference,1)==n && size(targets,1)==n,'rp:Endpoints','Endpoint counts must match or be scalar.');
local=([targets,repmat(g.rxPlaneZ,n,1)]-g.pRis.')*g.RRis;
height=local(:,3)-g.zTop;
path=rp_layer_path(local(:,1:2)-reference,[repmat(g.hOut,n,1),max(height,1e-9)], ...
    [c.optics.nGlass,c.optics.nAir],c.numerics.bisectionIterations);
o=struct('q',path.q,'opl',path.opl,'A',zeros(2,2,n),'B',zeros(2,2,n), ...
    'C',zeros(2,2,n),'slope',zeros(n,2),'valid',height>0);
T=g.RRis(1:2,1:2).'; z=g.RRis(1:2,3);
q=reshape(path.q.',2,1,n); outer=pagemtimes(q,pagetranspose(q));
kg=path.kz(:,1); ka=path.kz(:,2); identity=repmat(eye(2),1,1,n);
J=identity.*reshape(g.hOut./kg+height./ka,1,1,n)+ ...
    outer.*reshape(g.hOut./kg.^3+height./ka.^3,1,1,n);
o.A=pagemldivide(J,identity);
E=T-pagemtimes(q./reshape(ka,1,1,n),z.');
o.B=pagemtimes(o.A,E);
o.C=pagemtimes(pagetranspose(E),o.B);
o.slope=path.q*T+ka*z.';
end
