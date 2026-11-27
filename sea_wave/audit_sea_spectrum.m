function report = audit_sea_spectrum(outputFolder)
%AUDIT_SEA_SPECTRUM Check corrected quadrature against independent integration.
% Includes the former U=12, fmin=0.1 case and frequency-grid convergence.
    folder = fileparts(mfilename('fullpath'));
    if nargin == 0, outputFolder = fullfile(folder,'audit'); end
    if ~isfolder(outputFolder), mkdir(outputFolder); end
    sea = get_default_sea_params();
    p = get_default_RIS_params();
    winds = [3,5,9,12]; values = zeros(4,6);
    for j = 1:4
        p.windSpeed = winds(j);
        [wave,spectrum] = RIS_wave_components(p);
        variance = sum(abs(wave.coefficient).^2)/2;
        exact = bandVariance(p,sea.neumannConstant);
        fine = p; fine.frequencyBins = 400;
        denseWave = RIS_wave_components(fine);
        denseVariance = sum(abs(denseWave.coefficient).^2)/2;
        values(j,:) = [winds(j),4*sqrt(variance),spectrum.retainedFullSpectrumVariance, ...
            variance/exact-1,denseVariance/exact-1,sum(spectrum.frequencyWeightsHz)];
        assert(abs(variance/exact-1)<0.002,'Coarse spectral energy error exceeds 0.2%%.');
        assert(abs(denseVariance/exact-1)<0.00003,'Refined spectral quadrature did not converge.');
        assert(abs(sum(spectrum.frequencyWeightsHz)-(p.fMaxHz-p.fMinHz))<1e-12);
    end
    p.windSpeed = 12; p.fMinHz = 0.1;
    [wave,spectrum] = RIS_wave_components(p);
    exact = bandVariance(p,sea.neumannConstant);
    fixedError = sum(abs(wave.coefficient).^2)/(2*exact)-1;
    % Reconstruct the removed endpoint bug only for this numerical comparison.
    f = unique(wave.omega)/(2*pi);
    S = sea.neumannConstant*pi/4*(2*pi*f).^(-6).*exp(-2*p.g^2./(p.windSpeed*2*pi*f).^2);
    extraVariance = 2*pi*sum(S([1,end]).*spectrum.frequencyWeightsHz([1,end]));
    formerError = fixedError+extraVariance/exact;
    assert(abs(fixedError)<0.002 && formerError>0.05);
    cases = array2table(values,'VariableNames',{'windMps','HsM','retainedVariance', ...
        'defaultVarianceRelativeError','refinedVarianceRelativeError','sumFrequencyWeightsHz'});
    report = struct('cases',cases,'fMinHz',sea.f_min, ...
        'formerU12VarianceRelativeError',formerError,'correctedU12VarianceRelativeError',fixedError);
    writetable(cases,fullfile(outputFolder,'corrected_spectrum_audit.csv'));
    save(fullfile(outputFolder,'corrected_spectrum_audit.mat'),'report');
    disp(cases);
    fprintf('U=12, [0.1,4] Hz: former error %.6f%% -> corrected %.6f%%\n', ...
        100*formerError,100*fixedError);
end

function variance = bandVariance(p,C)
    variance = integral(@(w) C*pi/4*w.^(-6).*exp(-2*p.g^2./(p.windSpeed*w).^2), ...
        2*pi*p.fMinHz,2*pi*p.fMaxHz,'AbsTol',1e-13,'RelTol',1e-11);
end
