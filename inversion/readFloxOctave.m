function [wl, numUp, numDown, txt] = readFloxOctave(path, date)
    % path = fullfile('input', 'dataset inversion');
    % date = '230603';

    pathDown = fullfile(path, sprintf('Incoming_radiance_FULL_%s.CSV', date));
    pathUp = fullfile(path, sprintf('Reflected_radiance_FULL_%s.CSV', date));

    down = importdata(pathDown, ';');
    numDown = down.data;
    wl = numDown(:, 1);
    txt = down.colheaders;
    txt = cellfun(@(x) x(2:(end-1)), txt, 'UniformOutput', false);

    up = importdata(pathUp, ';');
    numUp = up.data;

    % [wl, spec, tts] = parseFlox(wl, numUp, numDown, txt, date, lon, lat);

end

