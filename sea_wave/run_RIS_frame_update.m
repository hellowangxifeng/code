function report = run_RIS_frame_update(evidence)
%RUN_RIS_FRAME_UPDATE Recompute and verify sea_wave and W_RIS for the 1 m RIS.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root,'sea_wave'),fullfile(root,'W_RIS'),fullfile(root,'fig3-7'));
    if nargin==0, evidence=fullfile(root,'研究执行','RIS_1m更新_20261105'); end
    if ~isfolder(evidence), mkdir(evidence); end
    diary(fullfile(evidence,'recompute.log')); guard=onCleanup(@() diary('off'));
    p=get_default_RIS_params();
    assert(isequal([p.Lx p.Ly],[1 1]),'Update entry expects the requested 1 m glass.');
    fprintf('=== Flotation and calm/regular-wave audit ===\n');
    report=audit_RIS_flotation(p);
    fprintf('%s\n',jsonencode(report.flotation));
    fprintf('=== Independent sea spectrum audit ===\n');
    spectrum=audit_sea_spectrum(fullfile(evidence,'spectrum'));
    report.spectrumCases=table2struct(spectrum.cases);
    fprintf('=== sea_wave: all 12 wind/direction cases ===\n');
    seaSuite=run_RIS_wind_comparison();
    report.seaCases=height(seaSuite.summary);
    report.seaSummary=table2struct(seaSuite.summary);
    report.seaGeometryPassed=verifySeaMotions(seaSuite.motions,p);
    clear seaSuite
    motionSummary=summarize_RIS_motion_existing();
    report.motionStatisticsCases=height(motionSummary);
    fprintf('=== Pure sea surface companion ===\n');
    renderPureSea(root);
    fprintf('=== W_RIS: production, 1e6 photons per panel ===\n');
    s=wris_run('production');
    wris_plot(s);
    report.wrisValidation=rmfield(s.validation,{'transparent','powers'});
    report.wrisVersion=s.config.version;
    report.wrisPhotonsPerPanel=s.config.mc.photons;
    report.wrisPanels=numel(s.panels);
    report.wrisGrid=s.config.grid.bins;
    report.completedAt=char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
    report.allRunsCompleted=true;
    save(fullfile(evidence,'update_audit.mat'),'report');
    writeJson(fullfile(evidence,'update_audit.json'),report);
    writeJson(fullfile(s.config.output,'acceptance_audit.json'), ...
        struct('version',s.config.version,'completedAt',report.completedAt, ...
        'flotationPassed',report.passed,'physicalAssertionsPassed',s.validation.physicalAssertionsPassed, ...
        'glassSizeM',[p.Lx p.Ly],'gridBins',s.config.grid.bins, ...
        'photonsPerPanel',s.config.mc.photons,'panels',numel(s.panels), ...
        'totalMassKg',report.flotation.totalMassKg,'auditFile',fullfile(evidence,'update_audit.json')));
    fprintf('COMPLETE: %s\n',fullfile(evidence,'update_audit.json'));
end

function passed=verifySeaMotions(motions,p)
    for j=1:numel(motions)
        o=motions{j};
        assert(isequal([o.parameters.Lx o.parameters.Ly],[p.Lx p.Ly]));
        assert(o.validity.numericalSolutionFinite);
        assert(o.validity.maxEquationResidual<1e-8);
        assert(o.validity.maxRotationOrthogonalityError<1e-11);
        sides=vecnorm(o.corners(2,:,:)-o.corners(1,:,:),2,2);
        assert(max(abs(sides-p.Lx),[],'all')<1e-11);
        sides=vecnorm(o.corners(3,:,:)-o.corners(2,:,:),2,2);
        assert(max(abs(sides-p.Ly),[],'all')<1e-11);
        assert(max(abs(o.boardCenter-o.bottomCenter-p.thickness*o.normal),[],'all')<1e-11);
    end
    passed=true;
end

function writeJson(file,value)
    fid=fopen(file,'w','n','UTF-8');
    assert(fid~=-1,'RIS_update:Write','Cannot open %s.',file);
    guard=onCleanup(@() fclose(fid));
    fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end

function renderPureSea(root)
% Keep the standalone script variables out of the update report workspace.
    run(fullfile(root,'sea_wave','main3_windseasurface_0326.m'));
end
