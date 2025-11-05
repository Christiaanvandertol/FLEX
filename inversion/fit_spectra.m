function results = fit_spectra(measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants, method)

if ~exist('lsqnonlin') %#ok<EXIST>
    pkg load optim
end

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
opt = optimset('MaxIter', 100, 'TolFun', stoptol);%, 'DiffMinChange', 1E-4); % works for float32 input
% 'Display', 'iter');

%% function minimization
f = @(params)COST_4SAIL_multiple(params, 1, measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants,0);

if any(tab.include)  % analogy of any(include == 1)
    tic
    [paramsout,Resnorm,FVAL,exitflag,output,~,J]= lsqnonlin(f, params0, lb, ub, opt); %#ok<ASGLU>
    toc
else % skip minimization and get resuls of RTMo_lite run with initial  parameters (param0)
    paramsout = params0;
end
er = f(paramsout);

%   s = [measurement.sigmarefl(:); tab.uncertainty(tab.include)];

xCov = inv(J.'*J)*Resnorm/(numel(FVAL)-numel(paramsout)); %#ok<MINV>
stdPar = sqrt(diag(full(xCov)));
%  stdPar = abs((inv(J.'*J)) * J.' * s);

%% best-fitting parameters
results = struct();
tab.value(tab.include) = paramsout;

results.parameters = demodify_parameters(tab.value, tab.variable);

%% best-fittiing spectra
f = @(params)COST_4SAIL_multiple(params, 0, measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants,1,stdPar,method);

[~, RSCOPE, L2C,FSCOPE] = f(paramsout);
J2 = numericalJacobian(f,paramsout);

varDiagnostic       = J2*xCov*J2';
stdDiagnostic       = sqrt(diag(varDiagnostic));
%results.rad = rad;
L2C.fAPAR_unc       = stdDiagnostic(1);
L2C.fAPARchl_unc    = stdDiagnostic(2);
L2C.FQE_unc         = sqrt(L2C.FQE_unc.^2 + stdDiagnostic(3).^2);
L2C.sigmaF_unc      = stdDiagnostic(4:end);


results.residual    = sqrt(er'*er);
results.Jacobian    = J;
results.L2C         = L2C;
results.FSCOPE      = FSCOPE;
results.RSCOPE      = RSCOPE;

% uncertainties of diagnostic outputs

% results.rmse = rmse;
% results.refl_mod = reflSAIL;
% results.sif = fluo.SIF;
% results.sif_norm = fluo.SIFnorm;
% results.soil_mod = soil.refl_in_meas;
% results.exitflag = exitflag;
end
