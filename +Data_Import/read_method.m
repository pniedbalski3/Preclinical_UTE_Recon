function [traj,Method_Params] = read_method(path)
%% Function to Read a Bruker Method File and return any parameters of interest
%Return a structure containing all the data that we care about

%% Make sure we have the path we need - no longer cd to path
if nargin == 0
    path=uigetdir('C:\','Select Folder in Which Method File is Located');
end

%% Get Method File/Files
[theo_method,meas_method] = Data_Import.identify_method_files(path);

%% Read Method File to get measured trajectories
locx = 0;
locy = 0;
locz = 0;
locend = 0;
%Read measured method file to get trajectories
fid=fopen(char(fullfile(path,meas_method)));
methodRead=textscan(fid,'%s','delimiter','\n');
methodRead=methodRead{1};
for index=1:size(methodRead,1)
    testStr=char(methodRead{index});
    if contains(testStr,'##$Method')
        Method=(testStr(11:end));
        Method_Params.Sequence = Method;
        Method_Params.SequenceName = Method;
    end
end
if contains(Method_Params.Sequence,'ute','IgnoreCase',true)
    Seq = 'Radial';
    Method_Params.Sequence = Seq;
elseif contains(Method,'spiral','IgnoreCase',true)
    Seq = 'Spiral';
    Method_Params.Sequence = Seq;
else
    error('Cannot Determine Sequence Type')
end

for index=1:size(methodRead,1)
    testStr=char(methodRead{index});
    if contains(testStr,'##$PVM_TrajKx')
        locx = index; %Index where X Trajectory begins 
        Length_X = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajKy')
        locy = index; %Index where Y Trajectory begins (and X ends)
        Length_Y = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajKz')
        locz = index; %Index where Z Trajectory begins (and Y ends)
        Length_Z = str2num(testStr(15:end));
    end
    if contains(testStr,'##$PVM_TrajBx')
        loc_end = index; %Index where Z Trajectory ends
    end
    if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
        AcqShift=str2num(testStr(13:end));
    end
    if contains(testStr,'##$PVM_SpatDimEnum=') %Dimension
        Dims = testStr(20:end);
    end
end
%Check length of Trajectories - if longer than 1, we have measured
%trajectories and should read them here. Otherwise move on
if Length_X > 1 
    %Measured Trajectories should work the same for both Spiral and Radial,
    %so only need minimal changes here.
    Method_Params.Traj_Type = 'measured';
    %Read out Trajectory Measurements
    trajx = [];
    trajy = [];
    trajz = [];
    for index = (locx+1):(locy-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        kx = str2num(readOutStr(1:end));
        trajx = [trajx kx];
    end
    for index = (locy+1):(locz-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        ky = str2num(readOutStr(1:end));
        trajy = [trajy ky];
    end
    for index = (locz+1):(loc_end-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        kz = str2num(readOutStr(1:end));
        trajz = [trajz kz];
    end
    %PJN - I hardcoded to measure 10 extra points at the end of a
    %trajectory in radial, so remove those here
    if strcmp(Seq,'Radial')
        NumExtra = 10;
    elseif strcmp(Seq,'Spiral')
        NumExtra = 0;
        if contains(Dims,'2D')
            trajx = trajx(1:(length(trajx)/2));
            trajy = trajy(((length(trajy)/2)+1):end);
        else
            trajx = trajx(1:(length(trajx)/3));
            trajy = trajy(((length(trajy)/3)+1):(2*length(trajy)/3));
            trajz = trajz((2*length(trajz)/3+1):end);
        end
    end
    trajx = trajx(1:(end-NumExtra));
    trajy = trajy(1:(end-NumExtra));
    
    trajx_shift_val = mean(trajx(1:AcqShift));
    trajy_shift_val = mean(trajy(1:AcqShift));
    
    trajx = trajx - trajx_shift_val; %These have a non-zero value at the beginning - get rid of that here
    trajy = trajy - trajy_shift_val;
    if contains(Dims,'3D')
        trajz = trajz(1:(end-NumExtra));
        trajz_shift_val = mean(trajz(1:AcqShift));
        trajz = trajz - trajz_shift_val;
    else
        trajz = zeros(size(trajx));
    end
else
    Method_Params.Traj_Type = 'theoretical';
    fid=fopen(char(theo_method));
    methodRead=textscan(fid,'%s','delimiter','\n');
    methodRead=methodRead{1};
    if strcmp(Seq,'Radial')  %Calculate Theoretical Trajectories for Radial
        ReadGrad = 0;
        PhaseGrad = 0;
        SliceGrad = 0;
        FlybackYN = nan;
        for index=1:size(methodRead,1)
            testStr=char(methodRead{index});
            if contains(testStr,'##$GradShape') %Theoretical Gradient Shape
                GradShapeStart = index;
                GradShapePts = str2num(testStr((strfind(testStr,'=')+1):end));
                for ind1 = index:size(methodRead,1)
                    testStr2 = char(methodRead{ind1});
                    if contains(testStr2,'##$')
                        GradShapeEnd = index;
                        break;
                    end
                end
            end
            if contains(testStr,'##$GradRes') %Gradient Shape Resolution
                GradRes=str2num(testStr(12:end));
                Method_Params.GradRes = GradRes;
            end
            if contains(testStr,'##$RampPoints') %Number of points for Ramp Compensation
                RampPoints = str2num(testStr(15:end));
            end
            if contains(testStr,'##$PVM_EffSWh=') %Bandwidth
                Bandwidth = str2num(testStr(15:end));
                Method_Params.Bandwidth = Bandwidth;
                Dwell = 1000/Bandwidth;
                Method_Params.Dwell = Dwell;
            end
            if contains(testStr,'##$FlybackYN') %Flyback Yes or No
                FlybackYN = testStr(14:end);
            end
            if contains(testStr,'##$FlybackPts') %Flyback Yes or No
                FlybackPts = str2num(testStr(15:end));
            end
            if contains(testStr,'##$ReadGrad')
                ReadGrad = str2num(testStr(13:end));
            end
            if contains(testStr,'##$PhaseGrad')
                PhaseGrad = str2num(testStr(14:end));
            end
            if contains(testStr,'##$SliceGrad')
                SliceGrad = str2num(testStr(14:end));
            end
            if contains(testStr,'##$PVM_GradCalConst')
                GradCalConst = str2num(testStr(21:end));
            end
            if contains(testStr,'##$PVM_Matrix=') %Matrix Size
                Matrix = str2num(char(methodRead{index+1}));
                %% Calculation is based on element 1
                NPts = Matrix(1);
            end
            if contains(testStr,'##$ExtraPoints')
                XtraPts = str2num(testStr(16:end));
            end
            if contains(testStr,'##$XtraPts')
                XtraPts = str2num(testStr(12:end));
            end
            if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
                AcqShift=str2num(testStr(13:end));
            end
        end
        GradShape = [];
        %PJN - Edit this to accomodate Flyback - need to add Flyback Pts to
        %method
        %Get the actual number of points
        if strcmp(FlybackYN,'No') || isnan(FlybackYN)
            Tot_Num_Pts = NPts+RampPoints+AcqShift+XtraPts;%NumPts + RampPoints;
            Method_Params.NPts = Tot_Num_Pts;
        else
            Tot_Num_Pts = FlybackPts;
            Method_Params.NPts = FlybackPts;
        end
        %Get the gradient shape
        for index = (GradShapeStart+1):(GradShapeEnd-1)
            readOutStr = methodRead{index};
            readOutStr = Tools.check_method_line_read(readOutStr);
            gradval = str2num(readOutStr(1:end));
            GradShape = [GradShape gradval];
        end
        %Integrate the gradient
        GradTime = 0:GradRes:(GradRes*(length(GradShape)-1));
        TrajShape = cumtrapz(GradTime,GradShape);
      
        SampleTime = 0:Dwell:(Dwell*(Tot_Num_Pts-1));
        %Get these in Physical Units if the correct, updated version is
        %used
        if ReadGrad ~= 0
            TrajShapex = GradCalConst * TrajShape * ReadGrad/100 /1000;
            TrajShapey = GradCalConst * TrajShape * PhaseGrad/100 /1000;
            TrajShapez = GradCalConst * TrajShape * SliceGrad/100 /1000;
            trajx = interp1(GradTime,TrajShapex,SampleTime);
            trajy = interp1(GradTime,TrajShapey,SampleTime);
            trajz = interp1(GradTime,TrajShapez,SampleTime);
        else    
            %Get the trajectory points at actual sampling points
            traj = interp1(GradTime,TrajShape,SampleTime);
            trajx = traj;
            trajy = traj;
            trajz = traj;
        end
    elseif strcmp(Seq,'Spiral')  %Calculate Theoretical Trajectories for Spiral
        for index=1:size(methodRead,1)
            testStr=char(methodRead{index});
            if contains(testStr,'##$Gxarray')
                GxStart = index;
                GxPts = str2num(testStr(12:end));
            end
            if contains(testStr,'##$Gyarray')
                GyStart = index;
                GyPts = str2num(testStr(12:end));
            end
            if contains(testStr,'##$Gzarray')
                GzStart = index;
                GzPts = str2num(testStr(12:end));
            end
            if contains(testStr,'##$Rx')
                RxStart = index;
                RxPts = str2num(testStr(9:end));
            end
            if contains(testStr,'##$GradRes') %Gradient Shape Resolution
                GradRes=str2num(testStr(12:end));
            end
            if contains(testStr,'##$PVM_EffSWh=') %Bandwidth
                Bandwidth = str2num(testStr(15:end));
                Method_Params.Bandwidth = Bandwidth;
                Dwell = 1000/Bandwidth;
                Method_Params.Dwell = Dwell;
            end
            if contains(testStr,'##$Num_Sample_Pts') %Number of sampling points
                NumPts=str2num(testStr(19:end));
            end
        end
        Gx = [];
        Gy = [];
        Gz = [];
        for index = (GxStart+1):(GyStart-1)
            readOutStr = methodRead{index};
            readOutStr = Tools.check_method_line_read(readOutStr);
            Gxval = str2num(readOutStr(1:end));
            Gx = [Gx Gxval];
        end
        for index = (GyStart+1):(GzStart-1)
            readOutStr = methodRead{index};
            readOutStr = Tools.check_method_line_read(readOutStr);
            Gyval = str2num(readOutStr(1:end));
            Gy = [Gy Gyval];
        end
        for index = (GzStart+1):(RxStart-1)
            readOutStr = methodRead{index};
            readOutStr = Tools.check_method_line_read(readOutStr);
            Gzval = str2num(readOutStr(1:end));
            Gz = [Gz Gzval];
        end
        Method_Params.GradShape_x = Gx;
        Method_Params.GradShape_y = Gy;
        Method_Params.GradShape_z = Gz;
        GradTime = 0:GradRes:(GradRes*(length(Gx)-1));
        TrajX = cumtrapz(GradTime,Gx);
        TrajY = cumtrapz(GradTime,Gy);
        TrajZ = cumtrapz(GradTime,Gz);
        SampleTime = 0:Dwell:(Dwell*(NumPts-1));
        trajx = interp1(GradTime,TrajX,SampleTime);
        trajy = interp1(GradTime,TrajY,SampleTime);
        trajz = interp1(GradTime,TrajZ,SampleTime);
    end
end
%Now, we have a single projection for each dimension whether we are measuring trajectories or not
%Get the rest of the parameters that we need from the theoretical method
%file
%Save the shapes in the method params structure
Method_Params.Base_Shape_x = trajx;
Method_Params.Base_Shape_y = trajy;
Method_Params.Base_Shape_z = trajz;

%% Read Parameters from Method File
fid=fopen(char(fullfile(path,theo_method)));
methodRead=textscan(fid,'%s','delimiter','\n');methodRead=methodRead{1};

Method_Params.ReadGrad = 0; %Line to make sure this gets set In versions before V3 of the code, this wasn't written to method file
for index=1:size(methodRead,1)
   	testStr=char(methodRead{index});
    if contains(testStr,'##OWNER=') %Scan Date, Time, Location
        DateTime = methodRead{index+1};
        ScanDate = DateTime(4:13);
        ScanTime = DateTime(15:22);
        Method_Params.ScanDate = ScanDate;
        Method_Params.ScanTime = ScanTime;
        FileLoc = methodRead{index+2};
        Method_Params.File = FileLoc;
    end
    if contains(testStr,'##$AcqShift') %Number of points for acquisition shift
        AcqShift=str2num(testStr(13:end));
        Method_Params.AcqShift = AcqShift;
    end
    if contains(testStr,'##$NPro') %Number of Projections
        NPro=str2num(testStr(9:end));
        Method_Params.NPro = NPro;
    end
    if contains(testStr,'##$NumTEs=') %NumTEs
        NumTEs=str2num(testStr(11:end));
        Method_Params.NumTEs = NumTEs;
    end
    if contains(testStr,'##$EchoTimes=') %Echo Times (just the indices - read out at the end)
        EchoTimesStart = index;
        for ind1 = (index+1):size(methodRead,1)
            testStr2 = char(methodRead{ind1});
            if contains(testStr2,'##$')
                EchoTimesEnd = index;
                break;
            end
        end
    end
    if contains(testStr,'##$SlabSelectYN=') %Are we slab selective?
        SlabYN = testStr(17:end);
        Method_Params.SlabYN = SlabYN;
    end
    if contains(testStr,'##$SlabThickness=') %Slab Thickness
        SlabThickness = str2num(testStr(18:end));
        Method_Params.SlabThickness = SlabThickness;
    end
    if contains(testStr,'##$ProjPerTrig=') %Projections Per Trigger (Useful for Xenon)
        ProjPerTrig = str2num(testStr(16:end));
        Method_Params.ProjPerTrig = ProjPerTrig;
    end
    if contains(testStr,'##$PVM_RepetitionTime=') %Repetition Time
        TR = str2num(testStr(23:end));
        Method_Params.TR = TR;
    end
    if contains(testStr,'##$ProUndersampling=') %Undersampling
        USamp = str2num(testStr(21:end));
        Method_Params.USamp = USamp;
    end
    if contains(testStr,'##$PVM_NAverages=') %Averages
        Averages = str2num(testStr(18:end));
        Method_Params.Averages = Averages;
    end
    if contains(testStr,'##$PVM_NRepetitions=') %Repetitions
        Repetitions = str2num(testStr(21:end));
        Method_Params.Repetitions = Repetitions;
    end
    if contains(testStr,'##$Num_Sample_Pts') %Number of sampling points
        NumPts=str2num(testStr(19:end));
        Method_Params.NPts = NumPts;
    end
    if contains(testStr,'##$RampPoints') %Number of points for Ramp Compensation
        RampPoints = str2num(testStr(15:end));
        Method_Params.RampPoints = RampPoints;
    end
    if contains(testStr,'##$PVM_EffSWh=') %Bandwidth
        Bandwidth = str2num(testStr(15:end));
        Method_Params.Bandwidth = Bandwidth;
        Dwell = 1000/Bandwidth;
        Method_Params.Dwell = Dwell;
    end
    if contains(testStr,'##$FlybackYN') %Flyback Yes or No
        FlybackYN = testStr(14:end);
        Method_Params.FlybackYN = FlybackYN;
    end
    if contains(testStr,'##$FlybackPts') %Flyback Yes or No
        FlybackPts = testStr(15:end);
        Method_Params.FlybackPts = str2num(FlybackPts);
    end
    if contains(testStr,'##$VFAYN') %Flyback Yes or No
        VFAYN = testStr(10:end);
        Method_Params.VFAYN = VFAYN;
    end
    if contains(testStr,'##$ExcPulse1Enum') %Pulse Shape
        PulseShape = testStr(18:end);
        Method_Params.PulseShape = PulseShape;
    end
    if contains(testStr,'##$PVM_Nucleus1Enum=') %Nucleus
        Nucleus = testStr(21:end);
        Method_Params.Nucleus = Nucleus;
    end
    if contains(testStr,'##$RefPowCh1=') %Reference Power
        RefPow = str2num(testStr(14:end));
        Method_Params.RefPow = RefPow;
    end
    if contains(testStr,'##$PVM_FrqWork=') %Working Frequency
        Freq = str2num(char(methodRead{index+1}));
        Freq = Freq(1);
        Method_Params.Frequency = Freq;
    end
    if contains(testStr,'##$PVM_Matrix=') %Matrix Size
        Matrix = str2num(char(methodRead{index+1}));
        Method_Params.MatrixSize = Matrix;
    end
    if contains(testStr,'##$PVM_EffSWh=') %Bandwidth
        Bandwidth = str2num(testStr(15:end));
        Method_Params.Bandwidth = Bandwidth;
    end
    if contains(testStr,'##$PVM_AcquisitionTime=') %Acquisition Time
        AcqTime = str2num(testStr(24:end));
        Method_Params.AcqTime = AcqTime;
    end
    if contains(testStr,'##$SpoilerAmp=') %Spoiler Amplitude
        SpoilerAmp = str2num(testStr(15:end));
        Method_Params.SpoilerAmp = SpoilerAmp;
    end
    if contains(testStr,'##$AmtSpoiling=') %Amount of spoiling
        Spoiling = str2num(testStr(16:end));
        Method_Params.Spoiling = Spoiling;
    end
    if contains(testStr,'##$SpoilDur=') %Spoiler Duration
        SpoilDur = str2num(testStr(13:end));
        Method_Params.SpoilDur = SpoilDur;
    end
    if contains(testStr,'##$RewinderYesNo=') %Rewinder
        Rewinder = testStr(18:end);
        Method_Params.Rewinder = Rewinder;
    end
    if contains(testStr,'##$RampTime=') %Ramp Time
        RampTime = str2num(testStr(13:end));
        Method_Params.RampTime = RampTime;
    end
    if contains(testStr,'##$DummyScans=') %Dummy Scans
        NumDummies = str2num(testStr(19:end));
        Method_Params.Dummies = NumDummies;
    end
    if contains(testStr,'##$PVM_TriggerModule=') %Trigger on or Off
        Trigger = testStr(22:end);
        Method_Params.Trigger = Trigger;
    end
    if contains(testStr,'##$ExcPulse1=') %Excitation Pulse Parameters
        PulseParams = testStr(14:end);
        PulseIndex = 1;
        while ~contains(char(methodRead{index+PulseIndex}),'##$')
            PulseParams = [PulseParams, char(methodRead{index+PulseIndex})];
            PulseIndex = PulseIndex+1;
        end
        Method_Params.PulseParams = PulseParams;
    end
    if contains(testStr,'##$ExcPul=') %Excitation Pulse Parameters
        PulseParams = testStr(14:end);
        PulseIndex = 1;
        while ~contains(char(methodRead{index+PulseIndex}),'##$')
            PulseParams = [PulseParams, char(methodRead{index+PulseIndex})];
            PulseIndex = PulseIndex+1;
        end
        Method_Params.PulseParams = PulseParams;
    end
    if contains(testStr,'##$PVM_SpatResol') %Resolution
        Resolution = str2num(methodRead{index+1});
        Method_Params.Resolution = Resolution;
    end
    if contains(testStr,'##$PVM_Fov=') %Field of View
        FOV = str2num(methodRead{index+1});
        Method_Params.FOV = FOV;
    end    
    if contains(testStr,'##$DiffusionYN=') %Diffusion Encoding YN
        DiffusionYN = testStr(16:end);
        Method_Params.DiffusionYN = DiffusionYN;
    end    
    if contains(testStr,'##$Nbvalue=') %Diffusion Encoding YN
        Nbvalue = str2num(testStr(12:end));
        Method_Params.Nbvalue = Nbvalue;
    end  
    if contains(testStr,'##$Bvalues=') %Diffusion Encoding YN
        bvalues = str2num(methodRead{index+1});
        Method_Params.bvalues = bvalues;
    end 
    if contains(testStr,'##$BDSG_Delta=') %Diffusion Encoding YN
        Delta = str2num(testStr(15:end));
        Method_Params.Delta = Delta;
    end 
    if contains(testStr,'##$BDSG_delta=') %Diffusion Encoding YN
        delta = str2num(testStr(15:end));
        Method_Params.delta = delta;
    end
    if contains(testStr,'##$GradRes') %Gradient Shape Resolution
        GradRes=str2num(testStr(12:end));
        Method_Params.GradRes = GradRes;
    end
    if contains(testStr,'##$XenonGradientsYesNo=') %Diffusion Encoding YN
        XeGradYN = testStr(24:end);
        Method_Params.XeGradYN = XeGradYN;
    end 
    %Params exclusively for Spiral
    if contains(testStr,'##$Rx')
        RxStart = index;
        RxPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Px')
        PxStart = index;
        PxPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Sx')
        SxStart = index;
        SxPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Ry')
        RyStart = index;
        RyPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Py')
        PyStart = index;
        PyPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Sy')
        SyStart = index;
        SyPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Rz')
        RzStart = index;
        RzPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Pz')
        PzStart = index;
        PzPts = str2num(testStr(7:end));
    end
    if contains(testStr,'##$Sz')
        SzStart = index;
        SzPts = str2num(testStr(7:end));
        for ind1 = (index+1):size(methodRead,1)
            testStr2 = char(methodRead{ind1});
            if contains(testStr2,'##$')
                SzEnd = ind1;
                break;
            end
        end
    end
    if contains(testStr,'##$Hubs=') %SNumber of Hubs
        Hubs = str2num(testStr(9:end));
        Method_Params.Hubs = Hubs;
    end
    if contains(testStr,'##$USampType=') %Undersampling type
        USampType = str2num(testStr(14:end));
        Method_Params.USampType = USampType;
    end
    if contains(testStr,'##$USampParams=') %Undersampling Parameters
        USampParams = str2num(methodRead{index+1});
        Method_Params.USampParams = USampParams;
    end
    if contains(testStr,'##$IAcqTime=') %Acquisition Time
        SpROTime = str2num(testStr(13:end));
        Method_Params.SpROTime = SpROTime;
    end
    if contains(testStr,'##$PVM_SPackArrNSlices=') %Number of slices
        NSlices = str2num(methodRead{index+1});
        Method_Params.NSlices = NSlices;
    end
    if contains(testStr,'##$PVM_SpatDimEnum=') %Number of Dimensions
        Dims = testStr(20:end);
        Method_Params.Dims = Dims;
    end
    if contains(testStr,'##$PVM_EchoTime=') %Echo Time (For Spiral
        TE = str2num(testStr(17:end));
        Method_Params.TE = TE;
    end
    if contains(testStr,'##$PVM_SliceThick=') %Slice Thickness
        SliceThick = str2num(testStr(18:end));
        Method_Params.SliceThick = SliceThick;
    end
    if contains(testStr,'##$NInterleaves=') %Slice Thickness
        Interleaves = str2num(testStr(17:end));
        Method_Params.Interleaves = Interleaves;
    end
    if contains(testStr,'##$ReadGrad')
        Method_Params.ReadGrad = str2num(testStr(13:end));
    end
    if contains(testStr,'##$PhaseGrad')
        Method_Params.PhaseGrad = str2num(testStr(14:end));
    end
    if contains(testStr,'##$SliceGrad')
        Method_Params.SliceGrad = str2num(testStr(14:end));
    end
    if contains(testStr,'##$ExtraPoints')
        Method_Params.ExtraPoints = str2num(testStr(16:end));
    end
    if contains(testStr,'##$XtraPts')
        Method_Params.ExtraPoints = str2num(testStr(12:end));
    end
    if contains(testStr,'##$PostPoints=')
        PostPoints = str2num(testStr(15:end));
        Method_Params.PostPoints = PostPoints;
    end
    if contains(testStr,'##$PVM_TxCoilAmpScaling1=')
        NCoil = str2num(testStr(26:end));
        Method_Params.NCoil = NCoil;
    end
end
%Make sure Various fields get set:
if ~isfield(Method_Params,'NCoil')
    Method_Params.NCoil = 1;
end
if ~isfield(Method_Params,'NumTEs')
    Method_Params.NumTEs = 1;
end
if ~isfield(Method_Params,'Nbvalue')
    Method_Params.Nbvalue = 1;
end

if strcmp(Seq,'Radial')
    TEs = [];
    for i = (EchoTimesStart+1):(EchoTimesEnd-1)
        TEhold = methodRead{i};
        TEhold = str2num(TEhold);
        TEs = [TEs TEhold];
    end
    Method_Params.TE = TEs;
    %% Rotate Projections according to golden means acquisition
    phi1 = 0.46557123;
    phi2 = 0.6823278;
    gs = 1;
    gr = 1;
    gp = 1;

    r = zeros(1,NPro);
    p = zeros(1,NPro);
    s = zeros(1,NPro);
    %Rotation code from Jinbang's UTE sequence
    for i = 0:(NPro-1)
        kz = (i*phi1-floor(i*phi1))*2-1;
        ts = kz*gs;
        alpha = (i*phi2-floor(i*phi2))*2*pi;
        tr = sqrt(1-kz*kz)*cos(alpha)*gr;
        tp = sqrt(1-kz*kz)*sin(alpha)*gp;
        r(i+1) = tr;
        p(i+1) = tp;
        s(i+1) = -ts;
    end

    trajx_rot = zeros(length(trajx),NPro);
    trajy_rot = trajx_rot;
    trajz_rot = trajx_rot;
    for i = 1:NPro
        trajx_rot(:,i) = r(i)*trajx';
        trajy_rot(:,i) = p(i)*trajy';
        trajz_rot(:,i) = s(i)*trajz';
    end

    traj = cat(3,trajx_rot,trajy_rot,trajz_rot);
    traj = permute(traj,[3,1,2]);
    %Trajectories are in physical Units at this point: Convert to Pixel
    %Units
    if Method_Params.ReadGrad == 0 && strcmp(Method_Params.Traj_Type, 'theoretical')
        radius = squeeze(sqrt(traj(1,:,:).^2+traj(2,:,:).^2+traj(3,:,:).^2));   
        maxradius = max(max(radius));
        traj = traj/maxradius/2;
        %traj(3,:,:) = traj(3,:,:) /(FOV(3)/FOV(1)); %Need to make the slice field of view behave itself
    else
        kFOV_desired = 1./Method_Params.Resolution;
        kMax_desired = kFOV_desired/2;
        max_k = max(kMax_desired);
        traj = traj/max_k/2;
        
        traj(3,:,:) = traj(3,:,:) * (FOV(3)/FOV(1)); %PJN - Add to get slice field of view correct... seems to have issues.
    end
elseif strcmp(Seq,'Spiral')
    %Read in rotation values
    Rx = [];
    Ry = [];
    Rz = [];
    Px = [];
    Py = [];
    Pz = [];
    Sx = [];
    Sy = [];
    Sz = [];
    for index = (RxStart+1):(PxStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Rxval = str2num(readOutStr(1:end));
        Rx = [Rx Rxval];
    end
    for index = (PxStart+1):(SxStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Pxval = str2num(readOutStr(1:end));
        Px = [Px Pxval];
    end
    for index = (SxStart+1):(RyStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Sxval = str2num(readOutStr(1:end));
        Sx = [Sx Sxval];
    end
    for index = (RyStart+1):(PyStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Ryval = str2num(readOutStr(1:end));
        Ry = [Ry Ryval];
    end
    for index = (PyStart+1):(SyStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Pyval = str2num(readOutStr(1:end));
        Py = [Py Pyval];
    end
    for index = (SyStart+1):(RzStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Syval = str2num(readOutStr(1:end));
        Sy = [Sy Syval];
    end
    for index = (RzStart+1):(PzStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Rzval = str2num(readOutStr(1:end));
        Rz = [Rz Rzval];
    end
    for index = (PzStart+1):(SzStart-1)
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Pzval = str2num(readOutStr(1:end));
        Pz = [Pz Pzval];
    end
    for index = (SzStart+1):(SzEnd-1) %For the spiral sequence, The AcqShift Line is the cutoff - was previously labeled EchoTimesEnd... maybe go back and change this at some point
        readOutStr = methodRead{index};
        readOutStr = Tools.check_method_line_read(readOutStr);
        Szval = str2num(readOutStr(1:end));
        Sz = [Sz Szval];
    end
    traj = zeros(3,NumPts,NPro);
%     Sx = Sx*FOV(3)/FOV(1);
%     Sy = Sy*FOV(3)/FOV(1);
%     Sz = Sz*FOV(3)/FOV(1);
    for i = 1:length(Rx)
        traj(1,:,i) = Rx(i)*trajx + Ry(i)*trajy + Rz(i)*trajz; 
        traj(2,:,i) = Px(i)*trajx + Py(i)*trajy + Pz(i)*trajz; 
        traj(3,:,i) = Sx(i)*trajx + Sy(i)*trajy + Sz(i)*trajz; 
    end
    if strcmp(Method_Params.Traj_Type, 'theoretical')
        traj = traj/max(traj(:))/2.2223;
        if contains(Dims,'3D')
            traj(3,:,:) = traj(3,:,:)*FOV(3)/FOV(1);
        end
    else
        kFOV_desired = 1./(Method_Params.Resolution);%/res_sphere);
        kMax_desired = kFOV_desired/2;
        max_k = max(kMax_desired);
        traj = traj/max_k/2;
        if contains(Dims,'3D')
            traj(3,:,:) = traj(3,:,:)*FOV(3)/FOV(1);
        end
    end
end
