function device = wris_device(wavelength)
%WRIS_DEVICE Shared LC electrode geometry and ideal air-calibrated limits.
arguments
    wavelength (1,1) double {mustBePositive,mustBeFinite}
end
device.electrodeWidthM = 0.7e-6;
device.electrodeGapM = 0.3e-6;
device.electrodePitchM = device.electrodeWidthM+device.electrodeGapM;
% Nyquist phase-gradient reference: ~15.426 deg at lambda0=532 nm.
ratio = wavelength/(2*device.electrodePitchM);
assert(ratio<1,'wris:ElectrodePitch','Require wavelength < twice the electrode pitch.');
device.maxAngleDeg = asind(ratio);
device.angleStepDeg = 0.021;
device.limitShape = 'box'; % independent air-calibrated X/Y commands
end
