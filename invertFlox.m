function Out = invertFlox(path_input, path_FLOX,numchar,formati)

%clear all %#ok<CLALL>
% for compiler comment paths
restoredefaultpath
addpath src/RTMs
addpath src/supporting
addpath src/fluxes
addpath src/IO
addpath src/inversion

%% 3. inputs
%path_input = fullfile('input');
%pathFlox = fullfile('input/dataset inversion');
%[wl, spec, angles] = readFlox(pathFlox);

tab             = readtable('flox_setup');
VarNames        = tab.Var1;
VarValues       = tab.Var2;
VarDescription  = tab.Var3;
VarUnits        = tab.Var4;

for k = 1:length(tab.Var1)
    FLOX.(VarNames{k}) = VarValues(k);
end

Efiles          = dir([path_FLOX '/*inc*Fluo*']);
Lfiles          = dir([path_FLOX '/*refl*Fluo*']);

if ~isempty(Efiles)
    Efilename   = [path_FLOX '/' Efiles(1).name];
else
    error(['no irradiance file found in ' path_FLOX])
end
if ~isempty(Lfiles)
    Lfilename    = [path_FLOX '/' Lfiles(1).name];
else
    error(['no upwelling radiance file found in ' path_FLOX])
end

if nargin<4
    formati     = '%1c%d%1c%d';
end
if nargin<3
    numchar = 2;
end

[wl, E, t]      = readFXBox(Efilename,numchar,formati);
[~, piL]        = readFXBox(Lfilename,numchar,formati);
[y,m,dom,hr,minute,sec] = datevec(t);

Doyt = t-datenum(['1-Jan' year(1)]);
Doy = floor(Doyt);
time = 24*(Doyt-Doy);




[sza,~,~,saa]   = calczenithangle(Doy,time,0,0,FLOX.lon,FLOX.lat);
angles.sza      = sza;
angles.vza      = FLOX.vza;
angles.psi      = FLOX.vaa - saa;

r = piL/E;
r(r<0) = NaN;
r(r>1) = NaN;

% spec = spec(:, 1);  % just to speed up
% angles.tts = angles.tts(1);
keyboard
%tab = read_input_sheet(pathFlox);
tab = read_input_sheet(path_input);

%% 4. fixed inputs
constants = define_constants();
spectral = define_bands();

pathPCflu = fullfile(pathFlox, 'PC_flu.csv');
PCflu = readtable(pathPCflu, "NumHeaderLines", 2);
pcf = table2array(PCflu);
pcf = pcf(:, 2:5);

pathFluspectPar = fullfile('input/fluspect_parameters/Optipar2021_ProspectPRO_CX.mat');
load(pathFluspectPar)  % optipar struct appears

atmfile = fullfile('input/radiationdata/FLEX-S3_std.atm');
atmo = load_atmo(atmfile, spectral.SCOPEspec);

% for atmo reading
meteo.Rin = 600;
meteo.Rli = 300;
meteo.Ta = 20;


%% 5. inversion

method = 'spline';  % M2020a name
if verLessThan('matlab', '9.8')
    method = 'splines';
end
specWlS     = interp1(wl, spec, spectral.wlS, method, NaN);

measurement.refl = specWlS;
results = fit_spectra(measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants);

tab.retrieved = results.parameters;
writetable(tab, 'inversion_results.csv')

%% check
figure
plot(spectral.wlS, measurement.refl)
hold on
plot(spectral.wlS, results.rad.refl)
% plot(spectral.wlS, rad.refl - measurement.refl)
set(gca, 'XScale', 'log')


%% 6. uncertainty

% propagation of uncertainty in refl to uncertainty in variables?
J = results.J;
meas_std_fit = std(specWlS, 0, 2);
meas_std_fit(isnan(meas_std_fit)) = 0;
stdPar = abs((inv(J.'*J)) * J.' * meas_std_fit);

% propagation of uncertainty in variables to uncertainty in reflectance
in_covar_mat = diag(tab.uncertainty(tab.include));  % we do not have and do not need covariances
out_covar_mat = J * in_covar_mat * J';
std_refl = diag(out_covar_mat);  % on diagonal var are located

