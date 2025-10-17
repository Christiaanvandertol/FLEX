function results = fit_spectra(measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants, method)


    %% here modification of parameters should occur
    tab = modify_tab_parameters(tab);
        
    %% save prior (x0) values
    tab.x0 = tab.value;
    
    %% initial parameter values, and boundary boxes
    iparams = tab.include;
    params0 = tab.value(iparams);
    lb = tab.lower(iparams);
    ub = tab.upper(iparams);
    
    stoptol = 1E-2;  % 
    %opt = optimset('MaxIter', 120, 'TolFun', stoptol, ...
    %               'DiffMinChange', 1E-4); % works for float32 input
    %                % 'Display', 'iter');
    opt = optimset('MaxIter', 60, 'TolFun', stoptol, ...
                   'DiffMinChange', 1E-4); % works for float32 input
                    % 'Display', 'iter');
    

    %% function minimization
    f = @(params)COST_4SAIL_multiple(params, measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants,0);

    if any(tab.include)  % analogy of any(include == 1)
        tic
        [paramsout,~,~,exitflag,output,~,J]= lsqnonlin(f, params0, lb, ub, opt); %#ok<ASGLU>
        toc
    else % skip minimization and get resuls of RTMo_lite run with initial  parameters (param0)
        paramsout = params0;
    end
    
    s = [measurement.sigmarefl(:); tab.uncertainty(tab.include)];
    stdPar = abs((inv(J.'*J)) * J.' * s);

    %% best-fitting parameters
    results = struct();
    tab.value(tab.include) = paramsout;

    results.parameters = demodify_parameters(tab.value, tab.variable);

    %% best-fittiing spectra
        f = @(params)COST_4SAIL_multiple(params, measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants,1,stdPar,method);

    [er, RSCOPE, L2C,FSCOPE] = f(paramsout);

    %results.rad = rad;
    results.residual    = sqrt(er'*er);
    results.Jacobian    = J;
    results.L2C         = L2C;
    results.FSCOPE      = FSCOPE;
    results.RSCOPE      = RSCOPE;
    
    % results.rmse = rmse;
    % results.refl_mod = reflSAIL;
    % results.sif = fluo.SIF;
    % results.sif_norm = fluo.SIFnorm;
    % results.soil_mod = soil.refl_in_meas;
    % results.exitflag = exitflag; 
end