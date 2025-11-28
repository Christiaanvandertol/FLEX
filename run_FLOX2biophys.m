% this is an example of looping through FLOX input data.
% the main code is 'FLOX2biophys

%direc = 'example input\';
direc = 'c:\users\tol\projects\FRM4FLUO\';

subdir = dir([direc 'JB*']);

for k = 1:length(subdir)
    subsubdir = dir([direc subdir(k).name '\2*']);
    path_specfit = [direc 'specfit\' subdir(k).name '_SpecFitOut\'];
    for j = 1:length(subsubdir)
        path_FLOX = [direc subdir(k).name '\' subsubdir(j).name '\'];
        %path_specfit = path_FLOX; % I placed the SIF in the same folder as the Reflectance
        L2valdata = FLOX2biophys(path_FLOX,path_specfit,1);
        filename = ['output/' subdir(k).name '_' subsubdir(j).name];
        writeL2C(filename,L2valdata)
        save( [direc subdir(k).name '\' subsubdir(j).name '\' subdir(k).name '_' subsubdir(j).name '_L2valdata.mat'],'L2valdata')
        save( ['output\' subdir(k).name '_' subsubdir(j).name '_L2valdata.mat'],'L2valdata')
    end
end
