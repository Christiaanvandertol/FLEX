function [out, refl, L2C,FSCOPE] = COST_4SAIL_multiple(p, minimize, measurement, tab, angles, ...
    spectral, optipar, pcf, atmo, meteo, constants, calcnetflux,stdPar,method)

% COST_4SAIL_my
% RETURNS
%   er                          difference between modeled and measured reflectance + prior weight

%% create input structs from table
tab.value(tab.include) = p;
tab.value = demodify_parameters(tab.value, tab.variable);

soilpar = table_to_struct(tab, 'soil');
canopy = table_to_struct(tab, 'canopy');
leafbio = table_to_struct(tab, 'leafbio');
wpcf = table_to_struct(tab, 'sif');


%% leaf reflectance - Fluspect
canopy.nlayers  = 60;
nl              = canopy.nlayers;

mly.nly      = 1;
mly.pLAI     = 1;  % the value does not matter
mly.totLAI   = 1;  % the value does not matter
mly.pCab     = leafbio.Cab;
mly.pCca     = leafbio.Cca;
mly.pCdm     = leafbio.Cdm;
mly.pCw      = leafbio.Cw;
mly.pCs      = leafbio.Cs;
mly.pN       = leafbio.N;

leafbio.V2Z = 0;
leafbio.Cbc = 0;
leafbio.Cp = 1;
leafbio.fqe = 0.01;

leafopt = fluspect_mSCOPE(mly,spectral,leafbio,optipar, nl);
%leafopt.refl(:, spectral.IwlT) = 0.01;
%leafopt.tran(:, spectral.IwlT) = 0.01;

%% soil reflectance - BSM
soilemp.SMC   = 25;        % empirical parameter (fixed) [soil moisture content]
soilemp.film  = 0.015;     % empirical parameter (fixed) [water film optical thickness]
% soilspec.wl  = optipar.wl;  % in optipar range
soilspec.GSV  = optipar.GSV;
soilspec.Kw   = optipar.Kw;
soilspec.nw   = optipar.nw;

soil.refl = BSM(soilpar, soilspec, soilemp);
%soil.refl(spectral.IwlT) = 0.06;

%% canopy reflectance factors - RTMo

canopy.x        = (-1/nl : -1/nl : -1)';         % a column vector
canopy.xl       = [0; canopy.x];                 % add top level
canopy.nlincl   = 13;
canopy.nlazi    = 36;
canopy.litab    = [ 5:10:75 81:2:89 ]';   % a column, never change the angles unless 'ladgen' is also adapted
canopy.lazitab  = ( 5:10:355 );           % a row
canopy.hot      = 0.05;
% canopy.hot  = canopy.leafwidth/canopy.hc;
canopy.lidf     = leafangles(canopy.LIDFa, canopy.LIDFb);

options.lite = 1;
options.calc_vert_profiles = 0;
options.calcnetflux = calcnetflux;
refl    = NaN*(ones(length(spectral.wlS),length(angles.tts)));
%SIF     =  zeros(length(spectral.wlS),length(angles.tts));
for k = 1:length(angles.tts)
    angles_i.tts = angles.tts(k);
    angles_i.tto = angles.tto(k);
    angles_i.psi = angles.psi(k);
    % radk   = models.RTMo_lite(soil, leafopt, canopy, angles_i);

    [rad,gap] = RTMo(spectral,atmo,soil,leafopt,canopy,angles_i,constants,meteo,options);
    
    % rad.rdd(:,k) = radk.rdd;
    % rad.rsd(:,k) = radk.rsd;
    % rad.rdo(:,k) = radk.rdo;
    % rad.rso(:,k) = radk.rso;
    % rad.refl(:, k) = radk.refl;
    refl(:, k) = rad.refl;

    if ~minimize
        integr              = 'layers';

        Ps = gap.Ps(1:nl);
        Ph = (1-Ps);
        canopy.Pnsun_Cab    = canopy.LAI*meanleaf(canopy,rad.Pnu_Cab,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha_Cab    = canopy.LAI*meanleaf(canopy,rad.Pnh_Cab,'layers',Ph); % net PAR Cab shaded leaves (photons)
        
        canopy.Pnsun        = canopy.LAI*meanleaf(canopy,rad.Pnu,integr,Ps); % net PAR Cab sunlit leaves (photons)
        canopy.Pnsha        = canopy.LAI*meanleaf(canopy,rad.Pnh,'layers',Ph); % net PAR Cab shaded leaves (photons)

        canopy.Pntot_Cab    = canopy.Pnsun_Cab+canopy.Pnsha_Cab; % net PAR Cab leaves (photons)
        canopy.Pntot        = canopy.Pnsun+canopy.Pnsha; % net PAR Cab leaves (photons)

        IwlPAR              = spectral.IwlPAR;
        %P                  = 0.001 * Sint((rad.Esun_(IwlPAR)+rad.Esun_(IwlPAR)),spectral.wlS(IwlPAR));
        %L2C.fAPARchl(k)     = canopy.Pntot_Cab./P; %#ok<*AGROW>
        L2C.fAPARchl(k)     = canopy.Pntot_Cab./rad.PAR; %#ok<*AGROW>
        L2C.fAPAR(k)        = canopy.Pntot./rad.PAR; %#ok<*AGROW>
        
        etau            = 1+0*rad.Pnu;
        etah            = 1+0*rad.Pnh;
        rad             = RTMf(constants,spectral,rad,soil,leafopt,canopy,gap,angles_i,etau,etah);

        ep              = constants.A*ephoton(spectral.wlF'*1E-9,constants);
        phi             = interp1(spectral.wlS,optipar.phi,spectral.wlF',method);
        EoutFrc_        = 1E-3*leafbio.fqe*ep.*(canopy.Pntot_Cab*phi); %1E-6: umol2mol, 1E3: nm-1 to um-1
        %EoutFrc     = 1E-3*Sint(EoutFrc_,spectral.wlF);
        %EoutFrc_(EoutFrc_<.01) = NaN;
        L2C.sigmaF(:,k)      = pi*rad.LoF_./EoutFrc_;

        if isfield(measurement,'sif')
            sifintm         = Sint(measurement.sif(:,k),spectral.wlF);
            s_sifintm       = Sint(measurement.sif_unc(:,k),spectral.wlF);
            sifint          = Sint(rad.LoF_,spectral.wlF);
            s_sifint        = 0;                        % SCOPE output uncertainty here.
            L2C.FQE(k)      = sifintm/sifint*leafbio.fqe;
           % L2C.FQE_unc(k)  = L2C.FQE(k).*sqrt( (s_sifintm./sifintm).^2 + (s_sifint./sifint).^2);
        end
    end

    %  %rad(k) = radk; %#ok<AGROW>
      %% canopy fluorescence from PCA, in W m-2 sr-1
           SIFi= pcf * cell2mat(struct2cell(wpcf));
           SIF_PCA(:,k) = interp1(640-399:850-399,SIFi,spectral.wlS,'linear',0);
    %%    rad.SIF(:,k) = SIF(640-399:850-399);
end
%% Biophysical data products FLEX ('L2C')
if ~minimize
    %measured iPAR
    if size(measurement.Ein,2)>1
        Ein             = measurement.Ein_all;
        fAPARchl        = interp1(angles.time,L2C.fAPARchl,(measurement.t_all-floor(measurement.t_all))*24);
        sigmaF        = interp1(angles.time,L2C.sigmaF',(measurement.t_all-floor(measurement.t_all))*24);
    else
        Ein             = measurement.Ein;
        fAPARchl        = L2C.fAPARchl;
        sigmaF          = L2C.sigmaF;
    end
    
    ep              = constants.A*ephoton(spectral.wlS(IwlPAR)*1E-9,constants);
    P               = 1E3*Sint(Ein(IwlPAR,:)./ep,spectral.wlS(IwlPAR));

    L2C.APARchl     = fAPARchl.*P;
    L2C.LCC         = leafbio.Cab;
    L2C.LCAR        = leafbio.Cca;
    L2C.LAI         = canopy.LAI;

    L2C.LAIunc      = stdPar(strcmp(tab.variable,'LAI'));
    L2C.LCCunc      = stdPar(strcmp(tab.variable,'Cab'));
    L2C.LCARunc     = stdPar(strcmp(tab.variable,'Cca'));
    ep              = constants.A*ephoton(spectral.wlF'*1E-9,constants);

    FSCOPE          = leafbio.fqe * phi*1E-3.*ep.*(sigmaF.*L2C.APARchl)';
    %stdDiagn        = J2*xCov*J2';
end

%% calculate the difference between measured and modeled data
er1 = (refl - measurement.refl - SIF_PCA./measurement.Ein)./measurement.sigmarefl;
er1 = er1(~isnan(er1));

%% add extra weight from prior information
prior.Apm = tab.x0(tab.include);
prior.Aps = tab.uncertainty(tab.include);
er2 = (p - prior.Apm) ./ prior.Aps;

%% total
er = [er1(:); 3E-2* er2];

if minimize
    out = er;
else
    out = [mean(L2C.fAPAR), nanmean(L2C.fAPARchl), nanmean(L2C.FQE), nanmean(L2C.sigmaF,2)' ]'; %#ok<NANMEAN>
end
end
