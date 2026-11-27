function [wave,spectrum] = RIS_wave_components(p,quadrature)
%RIS_WAVE_COMPONENTS Project Neumann/cos^2 sea using a local seeded RNG.
% Uses trapezoidal endpoint half-weights in frequency and direction.
% The optional 'trapezoid' argument is retained for callers such as W_RIS.
    arguments
        p (1,1) struct
        quadrature (1,:) char {mustBeMember(quadrature,{'trapezoid'})} = 'trapezoid'
    end
    f = logspace(log10(p.fMinHz),log10(p.fMaxHz),p.frequencyBins).';
    th = linspace(-pi/2,pi/2,p.directionBins);
    df = [(f(2)-f(1))/2;(f(3:end)-f(1:end-2))/2;(f(end)-f(end-1))/2];
    dw = repmat(th(2)-th(1),size(th));
    dw([1,end]) = dw([1,end])/2;
    [om,theta] = ndgrid(2*pi*f,th);
    sea = get_default_sea_params();
    S = sea.neumannConstant*pi/4*om.^(-6).*exp(-2*p.g^2./(p.windSpeed^2*om.^2));
    D = 2/pi*cos(theta).^2;
    amplitude = sqrt(2*S.*D.*(2*pi*df*dw));
    stream = RandStream('mt19937ar','Seed',p.seed);
    phase = 2*pi*rand(stream,size(amplitude));
    beta = theta+pi/2-deg2rad(p.windAngleDeg);
    k = om.^2/p.g;
    wave = struct('omega',om(:),'kx',reshape(k.*cos(beta),[],1), ...
        'ky',reshape(k.*sin(beta),[],1),'coefficient',amplitude(:).*exp(1i*phase(:)));
    a = 2*p.g^2/p.windSpeed^2;
    retained = gammainc(a/(2*pi*p.fMinHz)^2,2.5)-gammainc(a/(2*pi*p.fMaxHz)^2,2.5);
    spectrum = struct('name','project Neumann S(omega), cos^2 spreading', ...
        'quadrature',quadrature,'retainedFullSpectrumVariance',retained, ...
        'frequencyWeightsHz',df,'directionWeightsRad',dw,'frequencyBandHz',[p.fMinHz,p.fMaxHz]);
end
