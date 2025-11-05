direc = 'c:\Users\tol\projects\FRM4FLUO\';

subdir = dir([direc 'JB*']);

for k = 1:1%length(subdir)
    subsubdir = dir([direc subdir(k).name '\2*']);
    path_specfit = [direc 'specfit\' subdir(k).name '_SpecFitOut\'];
    for j = 1:1%length(subsubdir)
        path_FLOX = [direc subdir(k).name '\' subsubdir(j).name '\'];
        L2C_retrieval = FLOX2biophys(path_FLOX,path_specfit,1);
        save( [direc subdir(k).name '\' subsubdir(j).name '\' subdir(k).name '_' subsubdir(j).name '_retrieval.mat'],'L2C_retrieval')
    end
end