function [wl, data, time] = readFXBox(filename,numchar,formati,ID)

fileID      = fopen(filename);
T           = readtable(filename, 'TreatAsEmpty', '#N/D');
fclose(fileID);
x           = table2array(T);
if nargin<2
    numchar = 4;
end

if nargin<3
    formati = '%3c%d%c%d%c%f%c';
end

format = [ '%' num2str(numchar) 'c'];
for n = 1:size(x,2)
 %   format =[format '%2c%d%c%d%c%d%c']; %#ok<AGROW>
    format =[format formati]; %#ok<AGROW>
end

fid = fopen(filename);
line = fgetl(fid);
z = sscanf(line,format);
fclose('all');
%hr = z(7:8:end);
%minute = z(9:8:end);
%sec = z(11:8:end);


if nargin<3
    % FIRST OPTION
    %    hr = z(8-4+numchar:9:end);
    %    minute = z(10-4+numchar:9:end);
    %    sec = z(12-4+numchar:9:end);
    %    hr = z(7:8:end);
    %    minute = z(9:8:end);
    %    sec = z(11:8:end);

    % SECOND OPTION
    hr = z(8:9:end);
    minute = z(10:9:end);
    sec = z(12:9:end);

    % THIRD OPTION
else
    switch nargin
        case 3
            DDMMYY = z(4:4:end);
            HHMINSEC = z(6:4:end);
            dom = floor(DDMMYY/1E4);
            mon = floor((DDMMYY-dom*1E4)/1E2);
            year = 2000+DDMMYY-dom*1E4-mon*1E2;
            hr = floor(HHMINSEC/1E4);
            minute = floor((HHMINSEC-hr*1E4)/1E2);
            sec = HHMINSEC-hr*1E4-minute*1E2;

        case 4

            year = z(ID(1):14:end); % 7, 9, 11, 13, 15, 17
            mon  = z(ID(2):14:end);
            dom  = z(ID(3):14:end);
            hr      = z(ID(4):14:end);
            minute      = z(ID(5):14:end);
            sec      = z(ID(6):14:end);
    end



    %         ymd = z(8:7:end);
    %         dom = floor(ymd/1E4);
    %         mon = floor(ymd/1E2-dom*1E2);
    %         year = 2000+ymd-dom*1E4-mon*1E2;
    %         hrmins = z(10:7:end);
    %         hr = floor(hrmins/1E4);
    %         minute = floor(hrmins/1E2-hr*1E2);
    %         sec = hrmins-hr*1E4-minute*1E2;

end


%hr = z(7-4+numchar:8:end);
%minute = z(9-4+numchar:8:end);
%sec = z(11-4+numchar:8:end);


wl = x(:,1);
%data = x(:,2:Last+1);
data = x(:,2:end);

%time = datenum(0,0,0,hr,minute,sec); %#ok<*SAGROW>

time = datenum(year,mon,dom,hr,minute,sec); %#ok<*SAGROW>


if abs(length(time)- size(data,2))>0
    data = data(:,1:length(time));
end

% 
% if nargin<3
% %    hr = z(8-4+numchar:9:end);
% %    minute = z(10-4+numchar:9:end);
% %    sec = z(12-4+numchar:9:end);
% %    hr = z(7:8:end);
% %    minute = z(9:8:end);
% %    sec = z(11:8:end);
% 
%     hr = z(8:9:end);
%     minute = z(10:9:end);
%     sec = z(12:9:end);
% 
% 
% else
%     %keyboard
%     year = z(ID(1):14:end); % 7, 9, 11, 13, 15, 17
%     mon  = z(ID(2):14:end);
%     dom  = z(ID(3):14:end);
%     hr   = z(ID(4):14:end);
%     minute = z(ID(5):14:end);
%     sec  = z(ID(6):14:end);
% end
% 
% 
% %hr = z(7-4+numchar:8:end);
% %minute = z(9-4+numchar:8:end);
% %sec = z(11-4+numchar:8:end);
% 
% 
% wl = x(:,1);
% %data = x(:,2:Last+1);
% data = x(:,2:end);
% keyboard
% time = datenum(0,0,0,hr,minute,sec); %#ok<*SAGROW>
% 
% if abs(length(time)- size(data,2))>0
%     data = data(:,1:length(time));        
% end
% keyboard