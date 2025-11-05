function Out = FLOX2biophys(path_FLOXdata,path_SIFdata,option,pathSCOPE, path_settings)
% FLOX2biophys retrieves biophysical parameters from FLOX data
% authors: Christiaan van der Tol (c.vandertoL@utwente.nl) and Egor
% Prikaziuk (e.prikaziuk@utwente.nl)
% date: 28 August 2025
% The development of this function was funded by ESA-FLEX-DISC
%
% usage:
% Out = FLOX2SCOPE(path_input, path_FLOX,numchar,formati)
% this function retrieved biophysical properties from the spectra measured
% with the FLOX instrument. With modification it can be applied to other
% spectrometers.
% the retrieval is carried out with numerical inversion of the SCOPE mode.
% input:
%   path_FLOXdata: path to the location of the data of the FLOX
%   path_settings (optional): path to the location of the settings Default
%   is current folder.


%% notes for translators of the code to other languages:
% 1. this function uses rad2deg. This is simpy a multiplication of the variable by 180/pi;

%% 1. paths
%clear all %#ok<CLALL>
% for compiler comment paths

if ~exist('pathSCOPE','var')
    pathSCOPE = 'SCOPE'; % relative path to find SCOPE
end
if ~exist('path_settings','var')
    path_settings = ''; % settings are here
end

restoredefaultpath
addpath([pathSCOPE '/src/RTMs'])
addpath([pathSCOPE '/src/supporting'])
addpath([pathSCOPE '/src/fluxes'])
addpath([pathSCOPE '/src/IO'])
addpath('inversion')

%% 2. definitions
constants       = define_constants();
spectral        = define_bands();
% change the spectral bands compared to SCOPE (lower spectral resolution)
spectral.wlSori = spectral.wlS;
spectral.wlS    = (400:4:900)';
spectral.wlE    = (400:5:750)'; % this is necessary, compatibility with Fluspect
spectral.wlP    = spectral.wlS;
spectral.wlPAR  = spectral.wlS(spectral.wlS>=400 & spectral.wlS<=700);  % PAR range
spectral.IwlP   = 1:length(spectral.wlP);
spectral.IwlT   = length(spectral.wlP)-2:length(spectral.wlP);
spectral.IwlPAR = find(spectral.wlS>=400 & spectral.wlS<=700)';  % PAR range

spectral.wlT    = 898:900; % dummy
% the following part is Matlab specific

method = 'spline';  % M2020a name

%% 3. some FLOX specific settings
% this describes the header line of the FLOX data file.
% These are used in readFXBox.m, which could be made language specific/smarter.
formati         = '%1c%d%1c%d';
numchar         = 2;


%% 4. read the settings for the retrieval
tab             = readInputSheet(path_settings);

% read the FLOX setup
[A,B,~,~] = textread( [path_settings 'flox_setup.csv'], '%s %d %s %s' ,'delimiter' , ','); %#ok<DTXTRD>
VarValues       = double(B);
VarNames        = A;
%VarDesription   = C;
%VarUnits        = D;
for k = 1:length(VarNames)
    FLOX.(VarNames{k}) = VarValues(k);
end

%% 4. read the FLOX data
Efiles          = dir([path_FLOXdata '*inc*rad*Full*']);
Lfiles          = dir([path_FLOXdata '*refl*rad*Full*']);
ufiles          = dir([path_FLOXdata '*refl*unc*Full*']);
SIFfiles        = dir([path_SIFdata '*FLOX_SIF_allmeas*.txt']);
SIFuncfiles     = dir([path_SIFdata '*FLOX_SIF_uncertainty_allmeas*.txt']);


if isempty(Efiles)
    error(['no irradiance file found in ' path_FLOXdata])
end
if isempty(Lfiles)
    error(['no upwelling file found in ' path_FLOXdata])
end
if isempty(ufiles)
    warning(['no reflectance uncertainty file found in ' path_FLOXdata, ' ,using dummy instead'])
end
if isempty(SIFfiles)
    warning(['no SIF data file found in ' path_SIFdata, ' ,not calculating FQE'])
end

[E, t, piL,r_unc,SIF,SIF_unc,tiSIF] = deal([]);
for fileno = 1:length(Efiles)
    Efilename   = [path_FLOXdata '/' Efiles(fileno).name];
    Lfilename    = [path_FLOXdata '/' Lfiles(fileno).name];
    if ~isempty(SIFfiles)
        if length(SIFfiles)==length(Efiles) || fileno == 1
            SIFfilename         = [path_SIFdata  SIFfiles(fileno).name];
            SIFuncfilename      = [path_SIFdata  SIFuncfiles(fileno).name];
            [wlSIF, SIFi,tiSIFi]    = readFXbox(SIFfilename);
            [~, SIF_unci]    = readFXbox(SIFuncfilename);

            SIF_unci(isnan(SIF_unci)) = .2*SIFi(isnan(SIF_unci));
            kk = find(~isnan(mean(SIFi, 'omitnan')));
            tiSIFi = tiSIFi(kk);
            SIFi = SIFi(:,kk);
            SIF_unci = SIF_unci(:,kk);
        end
    end
    if ~isempty(ufiles)
        ufilename       = [path_FLOXdata '/' ufiles(fileno).name];
    end

    [wl, Ei, ti]        = readFXbox(Efilename,numchar,formati);
    [~, piLi]           = readFXbox(Lfilename,numchar,formati);
    if ~isempty(ufiles)
        [~, r_unci]        = readFXbox(ufilename,numchar,formati);
        r_unc   = [r_unc r_unci];
    end
    t       = [t; ti];
    E       = [E  Ei];
    piL     = [piL  piLi];
    if ~isempty(SIFfiles)
        if length(SIFfiles)==length(Efiles)
            SIF     = [SIF SIFi];
            SIF_unc = [SIF_unc SIF_unci];
            tiSIF   = [tiSIF tiSIFi];
        else
            SIF     = SIFi;
            SIF_unc = SIF_unci;
            tiSIF   = tiSIFi;
        end
    end
end

if ~isempty(SIFfiles)
    sif_m           = interp1(wlSIF,SIF,spectral.wlF,method,0);
    sif_unc         = interp1(wlSIF,SIF_unc,spectral.wlF,method,0);
    [tiSIF2,I]      = unique(tiSIF);
    sif             = sif_m(:,I)';
    sif_unc         = sif_unc(:,I)';
    sif             = interp1(tiSIF2,sif,t)';
    sif_unc         = interp1(tiSIF2,sif_unc,t)';
    calcFQE = 1;
else
    calcFQE = 0; % uncomment following lines for debugging only
end
r               = piL./E; % reflectance
r(r<0)          = NaN; % filtering, this is useful at at the edges of the spectrum or in low light conditions
r(r>1)          = NaN;
if isempty(r_unc)
    r_unc           = 0.01+ .05*r;        % this is the uncertainty of the FLOX reflectance. This is a dummy value for now!!
end

I       = isnan(r);
J       = find(sum(I)<50);
if length(J)<size(I,2)
    warning(['dataset contains ' num2str(size(I,2)-length(J)) ' poor quality spectra'])
end

if ~isempty(J)
    r       = r(:,J);
    t       = t(J);
    E       = E(:,J);
    piL     = piL(:,J);
    r_unc   = r_unc(:,J);

    % Unit conversion
    E       = E*1E3;            % from Wm-2nm-1 to Wm-2um-1
    piL     = piL*1E3;          %#ok<NASGU> % from Wm-2nm-1 to Wm-2um-1

    %% 5. calculate the angularity of the measurement setup.
    % This is needed in order to account for the BRDF

    y               = datevec(t);
    Doyt            = t-datenum(['1-Jan-' num2str(y(1))]); %#ok<DATNM> The decimal Julian calender date
    Doy             = floor(Doyt); % the Julian calander date
    time            = 24*(Doyt-Doy)+FLOX.timezone; % the time of the day in UTC
    [sza_rad,~,~,saa_rad]   = calczenithangle(Doy,time,0,0,FLOX.lon,FLOX.lat);
    angles.tts      = min(85,rad2deg(sza_rad));
    angles.tto      = single(repmat(single(FLOX.vza),length(angles.tts),1));
    angles.psi      = FLOX.vaa - rad2deg(single(saa_rad));

    %% 6. fixed inputs
    % fluorescence principle components. The algorithm does not retrieve
    % fluorescence, but it needs true reflectance, and thus needs to correct for fluorescence.
    %  If true reflectance is aready input, then specify this in the settings and the below will not be used.

    pathPCflu       = fullfile(path_settings, 'PC_flu.csv');
    %PCflu           = csvread(pathPCflu);
    PCflu           = dlmread(pathPCflu,',',1,0); %#ok<DLMRD>
    pcf             = PCflu(2:end, 2:5);

    pathFluspectPar = fullfile([pathSCOPE '/input/fluspect_parameters/Optipar2021_ProspectPRO_CX.mat']);
    load(pathFluspectPar)  %#ok<LOAD> % optipar struct appears
    optipar         = resampleOptipar(optipar,spectral.wlS); %#ok<NODEF>

    atmfile         = fullfile([pathSCOPE '/input/radiationdata/FLEX-S3_std.atm']);
    atmo            = load_atmo(atmfile, spectral.SCOPEspec);
    atmo.M          = interp1(spectral.wlSori, atmo.M, spectral.wlS);
    % for atmo reading
    % the ratio of Rin/Rli could influence the results if very long wavelengths
    % are included in the inversion (>2.5 um).
    % algorithm is not sensitive to these!
    meteo.Rin       = 600; % necessary to run SCOPE, arbitrary value
    meteo.Rli       = 300; % necessary to run SCOPE, about 0.5 of meteo.Rin
    meteo.Ta        = 20;  % necessary to run SCOPE, take any realstic value

    %% 7. inversion


    refl     = interp1(wl, r, spectral.wlS, method, NaN);
    Ein      = interp1(wl, E, spectral.wlS, method, NaN);
    refl_unc = interp1(wl, r_unc, spectral.wlS, method, NaN);

    % the option is an input, whether to do the retrieval for each spectrum
    % separately or for all at once. When all at once is chosen, it
    % limits to 10 per day, taking an interpolation.
    switch option
        case 0      % run for every spectrum separately
            uDoy = Doy;
            for k = 1:length(Doy)
                day(1).measurement(k).refl = refl(:,k);
                day(1).measurement(k).sigmarefl = refl_unc(:,k);
                day(1).measurement(k).Ein = Ein(:,k);
                day(1).measurement(k).t_all = t(k);
                day(1).angles(k).tts = angles.tts(k);
                day(1).angles(k).tto = angles.tto(k);
                day(1).angles(k).psi = angles.psi(k);
                day(1).angles(k).time = (t(k)-floor(t(k)))*24;
                if calcFQE
                    day(1).measurement(k).sif = sif(:,k);
                    day(1).measurement(k).sif_unc = sif_unc(:,k);
                end
            end
        case 1      % run one retrieval per day
            uDoy = unique(Doy);
            %allangles = angles;
            for d = 1:length(uDoy)          % loop over the days
                I = find(Doy==uDoy(d));
                day(d).measurement.Ein_all = Ein(:,I);
                day(d).measurement.t_all = t(I);
                day(d).measurement.tts_all = angles.tts(I);
                if length(I)>10     % if more than 10 measurements available on this day
                    [~,J] = sort(Doy(I));
                    Interval = (t(I(J(1))): (t(I(J(end)))-t(I(J(1))))/10 :t(I(J(end))))';
                    x = movmean(refl(:,I(J)),floor(length(I))/10,2);
                    day(d).measurement.refl = (interp1(t(I(J)), x', Interval))';
                    x = movmean(refl_unc(:,I(J)).^2,floor(length(I))/10,2);
                    day(d).measurement.sigmarefl = sqrt((interp1(t(I(J)), x', Interval))');
                    x = movmean(Ein(:,I(J)),floor(length(I))/10,2);
                    day(d).measurement.Ein = (interp1(t(I(J)), x', Interval))';
                    day(d).angles.tts = interp1(t(I(J)), angles.tts(I(J)), Interval);
                    day(d).angles.psi = interp1(t(I(J)), angles.psi(I(J)), Interval);
                    day(d).angles.tto = repmat(angles.tto(1),length(Interval),1); % this doesn't change
                    day(d).angles.time = (Interval-floor(Interval))*24;
                    if calcFQE
                        %x = movmean(sif(:,I(J)),floor(length(I))/10,2);
                        %day(d).measurement.sif = (interp1(t(I(J)), x', Interval))';
                        day(d).measurement.sif = sif(:,I);
                        %x = movmean(sif_unc(:,I(J)),floor(length(I))/10,2);
                        %day(d).measurement.sif_unc = interp1(t(I(J)), x', Interval))';
                        day(d).measurement.sif_unc = sif_unc(:,I);

                    end
                else                % if less than 11 measurements are available on this day
                    day(d).measurement.refl = refl(:,I);
                    day(d).measurement.sigmarefl = refl_unc(:,I);
                    day(d).measurement.Ein = Ein(:,I);
                    day(d).angles.tts = angles.tts(I);
                    day(d).angles.tto = angles.tto(I);
                    day(d).angles.psi = angles.psi(I);
                    day(d).angles.time = (t(I)-floor(t(I(1))))*24;
                    if calcFQE
                        day(d).measurement.sif = sif(:,I);
                        day(d).measurement.sif_unc = sif_unc(:,I);
                    end
                end
            end
    end

    for d = 1:length(day)
        measurement = day(d).measurement;
        angles = day(d).angles;
        for k = 1:length(day(d).angles)
            out = fit_spectra(measurement(k), tab, angles(k), ...
                spectral, optipar, pcf, atmo, meteo, constants, method);
            day(d).results(k).L2biophys = out.L2biophys; %#ok<*AGROW>
            day(d).results(k).FSCOPE = out.FSCOPE;
            day(d).results(k).RSCOPE = out.RSCOPE;
            day(d).results(k).L2biophys.sza = day(d).angles.tts;
            day(d).results(k).L2biophys.dec_hrs = day(d).angles.time;
            day(d).results(k).L2biophys.time_full = datestr(day(d).measurement(k).t_all); %#ok<DATST>
            day(d).results(k).wlR = spectral.wlP;
            day(d).results(k).wlF = spectral.wlF;
            x = day(d).angles.time;
            day(d).results(k).L2biophys.t = datestr(datenum(['1-Jan-' num2str(y(1))]) + uDoy(d) + x/24); %#ok<DATST,DATNM>
        end
        day(d).metadata=tab;
    end
    Out = day;
end
