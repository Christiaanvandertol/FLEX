function [wl, numUp, numDown, txt] = readFloxMat(path, date)
    % path = fullfile('input', 'dataset inversion');
    % date = '230603';

    pathDown = fullfile(path, sprintf('Incoming_radiance_FULL_%s.CSV', date));
    pathUp = fullfile(path, sprintf('Reflected_radiance_FULL_%s.CSV', date));
    
    down = readtable(pathDown);
    numDown = table2array(down);
    % wl = down.wl_f;
    wl = numDown(:, 1);  % for consistency with Octave
    txt = down.Properties.VariableNames;
   
    up = readtable(pathUp);
    numUp = table2array(up);

end
