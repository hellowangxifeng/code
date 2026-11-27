function field = bsf_phase_input(c, g, x, y, phiRis)
%BSF_PHASE_INPUT Ballistic scalar Gaussian phase under an explicit coherence assumption.
% Intensities are point samples here, not pixel averages. Scatter rays must be
% propagated as powers separately, never as sqrt(Is) with a common phase.
arguments
    c (1,1) struct
    g (1,1) struct
    x {mustBeNumeric, mustBeReal, mustBeFinite}
    y {mustBeNumeric, mustBeReal, mustBeFinite}
    phiRis {mustBeNumeric, mustBeReal, mustBeFinite} = 0
end
assert(isequal(size(x),size(y)), 'bsf:Shape', 'x and y must have matching sizes.');
assert(isscalar(phiRis)||isequal(size(phiRis),size(x)), 'bsf:Phase', 'RIS phase must be scalar or match x.');
q = [x(:),y(:),zeros(numel(x),1)]*g.RRis.'+g.pRis.';
r = (q-g.pTx.')*g.RTx;
z = r(:,3);
rho2 = sum(r(:,1:2).^2,2);
k = 2*pi*c.source.mediumRefractiveIndex/c.source.wavelength;
zRayleigh = pi*c.source.mediumRefractiveIndex*c.source.waist^2/c.source.wavelength;
inverseRadius = z./(z.^2+zRayleigh^2);
phase = reshape(k*z+k*rho2.*inverseRadius/2-atan(z/zRayleigh),size(x));
Ib = bsf_ballistic(c,g,x,y);
aperture = abs(x)<=c.ris.Lx/2 & abs(y)<=c.ris.Ly/2;
Ib = Ib.*aperture;
carrierGradient = k*(g.dTx.'*g.RRis(:,1:2));
field = struct('Ib',Ib,'phaseIncident',phase,'phaseRis',phiRis, ...
    'aperture',aperture, ...
    'Uincident',sqrt(Ib).*exp(1i*phase), ...
    'UafterRis',sqrt(Ib).*exp(1i*(phase+phiRis)), ...
    'carrierGradientRadPerM',carrierGradient, ...
    'carrierNyquistPitchM',pi./abs(carrierGradient), ...
    'gridPhaseSamplingCertified',false, ...
    'assumption','scalar coherent fundamental Gaussian ballistic field; incoherent scatter ensemble');
end
