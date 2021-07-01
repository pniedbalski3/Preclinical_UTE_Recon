function [theo_method,meas_method] = identify_method_files(path)

%Function to identify which method file contains measured trajectories
%using various clues from the method files.
%Input: path - path where data is stored

disp('Identifying Method Files')

%This should capture all method files present
methodfiles = dir(fullfile(path,'*ethod*'));

if length(methodfiles) > 2 %if 3 or more present, will have to guess about this
    disp('Too many method files present - Just guessing the correct method files')
elseif length(methodfiles) == 1 %If just 1 method file, use measured and theoretical
    theo_method = methodfiles.name;% fullfile(path,methodfiles.name);
    meas_method = methodfiles.name;% fullfile(path,methodfiles.name);
    return;
end

%Now do a bunch of comparisons of method files:
% 1) file size
% 2) date
% 3) Number of elements in ##$PVM_TrajKx
% 4) Value of ##$TrajUpToDate
% 5) Value of ##$TrajAdjMode

%I need the date of the raw data file to help identify the theoretical
%method file
%FID file might be named ser, rawdata.job0, or fid
rawfile = dir(fullfile(path,'rawdata.job0'));
if isempty(rawfile)
    rawfile = dir(fullfile(path,'fid'));
end
if isempty(rawfile)
    rawfile = dir(fullfile(path,'ser'));
end

rawfile_date = datenum(rawfile.date);

date_diff = zeros(1,length(methodfiles));
bytesize = zeros(1,length(methodfiles));
NumElements = zeros(1,length(methodfiles));
Trajup2Date = zeros(1,length(methodfiles));
TrajAdjMode = zeros(1,length(methodfiles));
MethodName = cell(1,length(methodfiles));
for i = 1:length(methodfiles)
    date_diff(i) = abs(datenum(methodfiles(i).date)-rawfile_date);
    bytesize(i) = methodfiles(i).bytes;
    fid=fopen(char(fullfile(path,methodfiles(i).name)));
    methodRead=textscan(fid,'%s','delimiter','\n');
    methodRead=methodRead{1};
    Length_X = nan;
    TrajUpToDate = nan;
    Mode = nan;
    Method = nan;
    for index=1:size(methodRead,1)
        testStr=char(methodRead{index});
        if contains(testStr,'##$Method')
            Method=(testStr(11:end));
        end
        if contains(testStr,'##$PVM_TrajKx')
            Length_X = str2num(testStr(15:end));
        end
        if contains(testStr,'##$PVM_TrajUpToDate')
            TrajUpToDate = testStr(21:end);
        end
        if contains(testStr,'##$PVM_TrajAdjMode')
            Mode = testStr(20:end);
        end
    end
    NumElements(i) = Length_X;
    MethodName{i} = Method;
    if contains(TrajUpToDate,'Yes')
        Trajup2Date(i) = 1;
    else
        Trajup2Date(i) = 0;
    end
    if contains(Mode,'Skip')
        TrajAdjMode(i) = 0;
    else
        TrajAdjMode(i) = 1;
    end
end
%Now, Compare dates - values that hint at theoretical go to 0. values that
%hint at measured go to 1
date_min = min(date_diff);
%Check if there are multiple with min value:
datemins = date_diff == date_min;
Date = ones(1,length(methodfiles)); %Theoretical should have values of 0
Date(datemins) = 0;

%Compare files size - usually measured should be bigger
byte_min = min(bytesize);
%Check if there are multiple with min value:
bytemins = bytesize==byte_min;
Bytes = ones(1,length(methodfiles)); %Theoretical should have values of 0
Bytes(bytemins) = 0;

%Compare number of points in the trajectory
numelmin = min(NumElements);
numelmins = numelmin == NumElements;
NumEl = ones(1,length(methodfiles)); %Theoretical should have values of 0
NumEl(numelmins) = 0;

%I also have the method name, which could be useful in ruling out a method
%file if more than 2 are provided... worry about that some other time.

Final = Date+Bytes+NumEl+Trajup2Date+TrajAdjMode;

% The smallest number has the most clues pointing to theoretical method.
[~,theoind] = min(Final);
% The largest number has the most clues pointing to measured method.
[~,measind] = max(Final);

theo_method = methodfiles(theoind).name;%fullfile(path,methodfiles(theoind).name);
meas_method = methodfiles(measind).name;

    
    
    