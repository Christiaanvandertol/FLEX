function [wl, spec, angles] = readFlox(path, wlName, specName)

    switch nargin 
        case 0
            path = fullfile('../../input/dataset inversion');
            wlName = 'wl_flox.csv';
            specName = 'Raveflox016m220809.csv';
        case 1
            wlName = 'wl_flox.csv';
            specName = 'Raveflox016m220809.csv';
    end
    
    wl = readtable(fullfile(path, wlName));
    wl = table2array(wl);
    wl = wl(5:end);  % I have wrong wl at the moment, quick solution

    spec = readtable(fullfile(path, specName), "NumHeaderLines", 4);
    spec = table2array(spec);

    [~, c] = size(spec);
    angles.tts = repelem(30, c);
    angles.tto = repelem(0, c);
    angles.psi = repelem(0, c);
end

