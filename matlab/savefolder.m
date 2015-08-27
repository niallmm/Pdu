function foldername = savefolder()

% this function creates part of a folder name with the date and time 

t = datetime;

foldername = datestr(t, 30);
% 
% year = num2str(t.Year);
% month = num2str(t.Month);
% day = num2str(t.Day);
% hour = num2str(t.Hour);
% min = num2str(t.Minute);

% foldername = [year '_' month '_' day '_' hour '_' min];