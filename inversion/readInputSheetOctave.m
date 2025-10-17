function tab = readInputSheetOctave(pathInput)
    tab = struct();
    % pathInput = fullfile('input', 'inversion_input_data.csv');
    data = importdata(pathInput);
    num = data.data;
    txt = data.textdata;
    
    [r, ~] = size(num);
    tab.tune = num(:, 1);
    tab.include = (tab.tune == 1);
    tab.variable = txt(2: (r + 1));
    tab.value = num(:, 3);
    tab.lower = num(:, 4);
    tab.upper = num(:, 5);
    tab.uncertainty = num(:, 6);
end