function wave = fx_wave_from_gauge(t,eta,opts)
%FX_WAVE_FROM_GAUGE Build RIS_float wave components from a measured gauge.
%   wave = FX_WAVE_FROM_GAUGE(t,eta) takes a uniformly sampled wave-elevation
%   record (m) and returns the deep-water components accepted by
%   RIS_float's p.wave: omega, kx, ky, coefficient.
%
%   The record is detrended and Fourier transformed. Positive-frequency bins
%   in [fMin,fMax] with amplitude above minAmplitude are retained; the wave is
%   assumed to propagate along +x. The synthesis reproduces the measured
%   elevation at x=0 to the FFT truncation error.
%
%   opts.fMin (0.15 Hz), opts.fMax (3 Hz), opts.minAmplitude (2e-4 m),
%   opts.tStart/opts.tEnd (empty = whole record).
arguments
    t (:,1) double {mustBeReal,mustBeFinite}
    eta (:,1) double {mustBeReal,mustBeFinite}
    opts.fMin (1,1) double {mustBePositive} = 0.15
    opts.fMax (1,1) double {mustBePositive} = 3
    opts.minAmplitude (1,1) double {mustBeNonnegative} = 2e-4
    opts.tStart (1,1) double = -Inf
    opts.tEnd (1,1) double = Inf
end
assert(numel(t)==numel(eta),'fx:Gauge','t and eta must have the same length.');
sel=t>=opts.tStart & t<=opts.tEnd;
tt=t(sel)-t(find(sel,1)); e=eta(sel)-mean(eta(sel));
fs=1/median(diff(tt)); N=numel(e); T=N/fs;
X=fft(e); f=(0:N-1).'/T; amp=2*abs(X)/N; ph=angle(X);
keep=f>opts.fMin & f<opts.fMax & amp>opts.minAmplitude;
wave=struct('omega',2*pi*f(keep),'kx',(2*pi*f(keep)).^2/9.81, ...
    'ky',zeros(nnz(keep),1),'coefficient',amp(keep).*exp(1i*ph(keep)));
wave.source='measured wave gauge, unidirectional +x, deep-water dispersion';
wave.fs=fs; wave.segment=[tt(1),tt(end)]; wave.components=nnz(keep);
wave.varianceM2=sum(abs(wave.coefficient).^2)/2;
wave.hsM=4*sqrt(wave.varianceM2);
end
