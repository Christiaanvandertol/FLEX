% this is an example of looping through FLOX input data.
% the main code is 'L2C_retrieval

direc = 'example input\';

subdir = dir([direc 'JB*']);

for k = 1:1%length(subdir)
    subsubdir = dir([direc subdir(k).name '\2*']);
    %path_specfit = [direc 'specfit\' subdir(k).name '_SpecFitOut\'];
    for j = 1:1%length(subsubdir)
        path_FLOX = [direc subdir(k).name '\' subsubdir(j).name '\'];
        path_specfit = path_FLOX; % I placed the SIF in the same folder as the Reflectance
        L2valdata = FLOX2biophys(path_FLOX,path_specfit,0);
        filename = [direc subdir(k).name '\' subsubdir(j).name '\' subdir(k).name '_' subsubdir(j).name '_retrieval.nc'];
        %writeL2C(filename,L2C_retrieval)
        save( ['output\' subdir(k).name '_' subsubdir(j).name '_L2valdata.mat'],'L2valdata')
    end
end
r