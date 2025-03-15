%% MAVENupstreamV7
%
% Analyzes MAVEN SWIA, SWEA, MAG and STATIC(L3) data to identify plasma
% regions and computes solar wind moments from SWIA distributions
%
% This version fits a convolution of the SWIA response function with a 
% 2-species anisotropic bi-Kappa distribution and a pickup ring distribution
%
% Features separate anisotropies for H+ and He++
%
% Note that naming conventions and some attributes differ from the source 
% MAVEN datasets
%
% Implementation log:
%        - Implement the SWIA response function in fine fits
%        - Crop domain for high response resolutiona cold fits 
%        - Force high resolution response function for very low-Kappa fits
%        - Use previous good fits as initial guess (could use the best fit specifically)
%        - Ingest a list of manual (override) fit parameters (beam direction implemented, might also need resp. funcition selection etc.)
%        - Calculate uncertainties on VHx,VHy,VHz and VAx,VAy,VAz as well
%          as THperp, THpara, TAperp, TApara
%        - Calculate total TH+ temperature (treating core and beam as one distribution)
%        - Attenuator-dependent response matrix selection
%        - Constraining the beam to SWIA fine H+ velocity space (v7.1)
%        - Adaptive selection of response matricies with wide sampling in azimuth for high elevations ('v7.2')
%
% Internal version: 7.2
%
% Author: Robin Ramstad (LASP/University of Colorado Boulder) 2023-2025

disp(' % MAVENupstreamV7');

clear;

%Search paths
saveLoc = '.\UpstreamFilesV7\';
savePlotsLoc = '.\UpstreamProcPlotsV7\';

%Settings
fcorrSWIA = 1; %Correction factor for SWIA diff. fluxes (use 1/0.7 with SWIA V1 data, currently orbits 1-16000)

Orbits = 1:22318; %2014 to Oct 2024
doplot = 1; % =1: Plot individual fits

%Physical constants
 Rm = 3390; %[km]
amu = 1.6605E-27; %[kg]
 mH = 1.008*amu;
 mA = 4.002*amu;
 mO = 16*amu;
mO2 = 32*amu;
 qH = 1; %H+ charge [e]
 qA = 2; %He2+ charge [e]
 me = 5.486E-4*amu; %Electron mass
 eV = 1.6027E-19;% [J/eV]
  k = 1.3806E-23;% [J/K] 
  q = 1.602176E-19; %Elementary charge [C] 
mu0 = 4*pi*1E-7; %Magnetic permeability in vacuum [H/m]
  
%Load list and times of MAVEN orbits
load('MAVENorbits.mat'); 

%If variable "Orbits" not set, search for unprocessed orbits
if ~exist('Orbits','var')
    disp('  %% Variable "Orbits" not set, searching for unprocessed orbits');
    Orbits = zeros(1,0);
    for orb = 1:max(OrbNum)
        if exist([SWIAL2matLoc,'mvn_swia_coarse_mso_orb',sprintf('%05d',orb),'.mat'],'file') && ~exist([saveLoc,'mvn_swi_sw_',sprintf('%05d',orb),'_v07_r01.mat'],'file')
            Orbits = cat(2,Orbits,orb);
        end
    end
    disp(['   %% Found ',num2str(length(Orbits)),' unprocessed orbits, processing']);
end

if max(Orbits) > max(OrbNum)
    error('  %% Maximum orbit number to convert is larger than orbit file, running  mergeOrbitFilesV2.m');
end

%Load SWIA response functions

    %Low resolution samples
    SWIArc{1} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDphInfDth20k1_n4.mat'); %Attenuator closed
    SWIAro{1} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDphInfDth20k1_n4.mat'); %Attenuator open

    %Mid-resolution samples
    SWIArc{2} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDph4Dth5k04_n12.mat'); %Attenuator closed
    SWIAro{2} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDph4Dth5k04_n12.mat'); %Attenuator open
    
    %High-resolution samples
    SWIArc{3} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDph2Dth2k02_n100v2.mat'); %Attenuator closed
    SWIAro{3} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDph2Dth2k02_n100v2.mat'); %Attenuator open
    
    %High-resolution samples w. wide range (for high elevations)
    SWIArc{4} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDph2Dth2k02_n200.mat'); %Attenuator closed
    SWIAro{4} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDph2Dth2k02_n400.mat'); %Attenuator open
    
    %Super-high-resolution samples (only for extremely low-kappa dists)
    SWIArc{5} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDph1Dth1k01_n400.mat'); %Attenuator closed
    SWIAro{5} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDph1Dth1k01_n400.mat'); %Attenuator open
    
    %Super-high-resolution samples w. wide range (only for extremely low-kappa dists at high elevations)
    SWIArc{6} = load('.\SWIAresponse\SWIAsims\SWIArespAttClosedDph1Dth1k01_n800.mat'); %Attenuator closed
    SWIAro{6} = load('.\SWIAresponse\SWIAsims\SWIArespAttOpenDph1Dth1k01_n800.mat'); %Attenuator open

%Suppress warnings from bad fits
warning('off','all');

indOrbs = find(ismember(OrbNum,Orbits));

% Load manually identified SW inclusion intervals
fid = fopen('SWincludeIntervals.txt'); kk = 0;
timeManInclSW = zeros(0,2);
str = fgetl(fid); %Run once to pass the header line
while str ~= -1
    kk = kk+1;
    str = fgetl(fid);
    if length(str) == 27
        timeManInclSW = cat(1,timeManInclSW,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
end
disp([' %% Found ',num2str(size(timeManInclSW,1)),' manually identified SW inclusion intervals in SWincludeIntervals.txt']);
fclose(fid);

% Load manually identified SW exclusion intervals
fid = fopen('SWexcludeIntervals.txt'); kk = 0;
timeManExclSW = zeros(0,2);
str = fgetl(fid); %Run once to pass the header line
while str ~= -1
    kk = kk+1;
    str = fgetl(fid);
    if length(str) == 27
        timeManExclSW = cat(1,timeManExclSW,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
end
disp([' %% Found ',num2str(size(timeManExclSW,1)),' manually identified SW exclusion intervals in SWexcludeIntervals.txt']);
fclose(fid);

%% Load manually set SW fitting model parameters intervals
fid = fopen('SWmodelParams.txt'); kk = 0;
timeManSWbeamSlow = zeros(0,2);
timeManSWbeamFast = zeros(0,2);
timeManSWbeamBneg = zeros(0,2);
timeManSWbeamBpos = zeros(0,2);
timeManDVA = zeros(0,2);
timeManTHcore = zeros(0,2);
timeManTATH = zeros(0,2);
timeManAlphaH = zeros(0,2);
timeManAlphaA = zeros(0,2);
timeManElim = zeros(0,2);
timeManAzLim = zeros(0,2);
timeManElLim = zeros(0,2);
timeManSWuseindr = zeros(0,2);
manDVA    = zeros(0,3);
manTHcore = zeros(0,3);
manTATH   = zeros(0,3);
manAlphaH = zeros(0,3);
manAlphaA = zeros(0,3);
manElim = zeros(0,2);
manAzLim = zeros(0,2);
manElLim = zeros(0,2);
manSWuseindr = zeros(0,1);
str = fgetl(fid); %Run once to pass the header line
while str ~= -1
    kk = kk+1;
    str = fgetl(fid);
    if ~strcmp(str(1),'%') && (~isempty(findstr(str,'beam s')) || ~isempty(findstr(str,'beam slow')))
        timeManSWbeamSlow = cat(1,timeManSWbeamSlow,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
    if ~strcmp(str(1),'%') && (~isempty(findstr(str,'beam f')) || ~isempty(findstr(str,'beam fast')))
        timeManSWbeamFast = cat(1,timeManSWbeamFast,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'beam -B'))
        timeManSWbeamBneg = cat(1,timeManSWbeamBneg,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'beam +B'))
        timeManSWbeamBpos = cat(1,timeManSWbeamBpos,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'THcore ')) 
        timeManTHcore = cat(1,timeManTHcore,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'THcore ')+7)))
            indCommas = findstr(str,',');
            temp = str2num(str((findstr(str,'THcore ')+7):(indCommas(find(sort([indCommas findstr(str,'THcore ')]) == findstr(str,'THcore ')))-1)));
            manTHcore = cat(1,manTHcore,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'TA/TH ')) 
        timeManTATH = cat(1,timeManTATH,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'TA/TH ')+6)))
            indCommas = findstr(str,',');
            temp = str2num(str((findstr(str,'TA/TH ')+6):(indCommas(find(sort([indCommas findstr(str,'TA/TH ')]) == findstr(str,'TA/TH ')))-1)));
            manTATH = cat(1,manTATH,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'alphaH ')) 
        timeManAlphaH = cat(1,timeManAlphaH,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'alphaH ')+7)))
            temp = [str2num(str(findstr(str,'alphaH ')+7)),str2num(str(findstr(str,'alphaH ')+9)),str2num(str(findstr(str,'alphaH ')+11))];
            manAlphaH = cat(1,manAlphaH,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'alphaA ')) 
        timeManAlphaA = cat(1,timeManAlphaA,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'alphaA ')+7)))
            indCommas = findstr(str,',');
            temp = str2num(str((findstr(str,'alphaA ')+7):(indCommas(find(sort([indCommas findstr(str,'alphaA ')]) == findstr(str,'alphaA ')))-1)));
            manAlphaA = cat(1,manAlphaA,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'Elim ')) 
        timeManElim = cat(1,timeManElim,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'Elim ')+5)))
            indCommas = findstr(str,',');
            temp = str2num(str((findstr(str,'Elim ')+5):(indCommas(find(sort([indCommas findstr(str,'Elim ')]) == findstr(str,'Elim ')))-1)));
            manElim = cat(1,manElim,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'AzLim ')) 
        timeManAzLim = cat(1,timeManAzLim,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'AzLim ')+6)))
            indCommas = findstr(str,',');
            temp = str2num(str((findstr(str,'AzLim ')+6):(indCommas(find(sort([indCommas findstr(str,'AzLim ')]) == findstr(str,'AzLim ')))-1)));
            manAzLim = cat(1,manAzLim,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'ElLim ')) 
        timeManElLim = cat(1,timeManElLim,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        indCommas = findstr(str,',');
        temp = str2num(str((findstr(str,'ElLim ')+6):(indCommas(find(sort([indCommas findstr(str,'ElLim ')]) == findstr(str,'ElLim ')))-1)));
        if ~isempty(temp)manElLim = cat(1,manElLim,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'dVA ')) 
        timeManDVA = cat(1,timeManDVA,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'dVA ')+(4:7))))
            temp = [str2num(str(findstr(str,'dVA ')+(4:7))),str2num(str(findstr(str,'dVA ')+(9:12))),str2num(str(findstr(str,'dVA ')+(14:17)))];
            manDVA = cat(1,manDVA,temp);
        end
    end
    if ~strcmp(str(1),'%') && ~isempty(findstr(str,'indr ')) 
        timeManSWuseindr = cat(1,timeManSWuseindr,[datenum(str(1:13),'yyyymmddTHHMM') datenum(str(15:27),'yyyymmddTHHMM')]);
        if ~isempty(str2num(str(findstr(str,'indr ')+5)))
            manSWuseindr = cat(1,manSWuseindr,str2num(str(findstr(str,'indr ')+5)));
        elseif length(str) > findstr(str,'indr ')+6
            if strcmp(str(findstr(str,'indr ')+(5:7)),'max')
                manSWuseindr = cat(1,manSWuseindr,length(SWIArc));
            end
        end
    end

end
if size(timeManSWuseindr,1) ~= size(manSWuseindr,1)
    error('Number of response selection intervals not equal to number of response selections (debuggin'' time!)');
end
disp([' %% Found ',num2str(size(timeManSWbeamSlow,1)),' manually identified slow-beam intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManSWbeamFast,1)),' manually identified fast-beam intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManSWbeamBneg,1)),' manually identified -B beam intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManSWbeamBpos,1)),' manually identified +B beam intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManDVA,1)),' (',num2str(size(manDVA,1)),' settings) manually set dVA constraint intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManTHcore,1)),' (',num2str(size(manTHcore,1)),' settings) manually set THcore constraint intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManTATH,1)),' (',num2str(size(manTATH,1)),' settings) manually set TA/TH constraint intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManAlphaH,1)),' (',num2str(size(manAlphaH,1)),' settings) manually set alpha_H constraint intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManAlphaA,1)),' (',num2str(size(manAlphaA,1)),' settings) manually set alpha_A constraint intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManElim,1)),' (',num2str(size(manElim,1)),' settings) manually set Elim intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManAzLim,1)),' (',num2str(size(manAzLim,1)),' settings) manually set Azimuth intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManElLim,1)),' (',num2str(size(manElLim,1)),' settings) manually set Elevation intervals in SWmodelParams.txt']);
disp([' %% Found ',num2str(size(timeManSWuseindr,1)),' (',num2str(size(manSWuseindr,1)),' settings) manually set response matricies intervals in SWmodelParams.txt']);
fclose(fid);

for indNum = 1:length(indOrbs) 
    tstartOrbit = tic;
    
    indOrb = indOrbs(indNum);
    
    disp([' %% Processing orbit #',num2str(OrbNum(indOrb)),' at ',datestr(now,'mmm dd - HH:MM:SS')]);
    
    %Load files
    clearvars SWIAc SWIAf STATIC SWEA MAG data;
    
    try SWIAc = load([SWIAL2matLoc,'mvn_swia_coarse_mso_orb',sprintf('%05d',OrbNum(indOrb)),'.mat']);
        SWIAc.j = fcorrSWIA.*SWIAc.j; %Apply SWIA correction factor
    catch
        warning(['No SWIA coarse file found for orbit ',sprintf('%05d',OrbNum(indOrb)),', skipping.']);
        continue;
    end
    
    try SWIAf = load([SWIAL2matLoc,'mvn_swia_fine_mso_orb',sprintf('%05d',OrbNum(indOrb)),'.mat']); SWIAfFileFound = 1;
        SWIAf.j = fcorrSWIA.*SWIAf.j; %Apply SWIA correction factor
        SWIAf.C = fcorrSWIA.*SWIAf.C; %Apply SWIA correction factor
    catch
        SWIAf.time = NaN(1,0);
        SWIAf.E = NaN(48,0);
        SWIAf.j = NaN(12,10,48,0);
        SWIAf.c = NaN(12,10,48,0);
        SWIAf.C = NaN(12,10,48,0);
        SWIAf.Mswia2msoAll = NaN(3,3,0);
        SWIAf.Mmso2swiaAll = NaN(3,3,0);
        SWIAf.ph = NaN(12,48,0);
        SWIAf.th = NaN(10,1);
        SWIAf.atten_state = NaN(1,0);
        SWIAf.xSC = NaN(1,0);
        SWIAf.ySC = NaN(1,0);
        SWIAf.zSC = NaN(1,0);
        warning(['No SWIA fine file found for orbit ',sprintf('%05d',OrbNum(indOrb)),', analyzing anyway.']);
        %continue;

        SWIAfFileFound = 0;
    end
    
    try SWEA  = load([SWEAL2matLoc,'mvn_swea_mso_orb',num2str(OrbNum(indOrb)),'.mat']); SWEAfileFound = 1;
    catch
        warning('Could not load SWEA file');
        SWEA.time = zeros(0,1);
        SWEA.j = zeros(6,16,64,0);
        SWEA.E = zeros(6,16,64,0);
        SWEA.xmso = zeros(6,16,64,0);
        SWEA.ymso = zeros(6,16,64,0);
        SWEA.zmso = zeros(6,16,64,0);
        SWEA.elev = NaN(6,64);
        SWEA.azim = NaN(16,1);
        SWEA.sc_pot = 0;
        
        SWEAfileFound = 0;
    end
    
    try STA  = load([STATICL3matLoc,'mvn_sta_mso_orb',sprintf('%05d',OrbNum(indOrb)),'_moments.mat']); STATICfileFound = 1;
    catch
        STA.time = zeros(0,1);
        STA.NH = zeros(0,1);
        STA.NO = zeros(0,1);
        STA.NO2 = zeros(0,1);
        warning('Could not load STATIC file');
        %STA.time = zeros(1,0);
        STATICfileFound = 0;
    end

    try SCpot = load(['D:\Data\MAVEN\SWEA\L3mat\scpotOrb\mvn_swe_l3_scpot_orb',sprintf('%05d',OrbNum(indOrb)),'.mat']); SCpotFileFound = 1;
        indRem = isnan(SCpot.pot);
        SCpot.time(indRem) = [];
        SCpot.pot(indRem) = [];
        SCpot.method(indRem) = [];
    catch
        SCpot.time = zeros(0,1);
        SCpot.pot = zeros(0,1);
        SCpot.method = zeros(0,1);
        SCpotFileFound = 0;
    end

    try MAG = load([MAGL2matLoc,'mvn_mag32Hz_mso_orb',sprintf('%05d',OrbNum(indOrb)),'.mat']);
    catch
        MAG.time = 0;
        MAG.Bx = NaN;
        MAG.By = NaN;
        MAG.Bz = NaN;
        warning('Could not load MAG file');
    end
    
    %% SWIA distributions and moments
    %Coarse
    SWIAc.dE = single(NaN(size(SWIAc.E)));
    for jj = 2:size(SWIAc.E,1)-1
        SWIAc.dE(jj,:) = (abs(SWIAc.E(jj-1,:)-SWIAc.E(jj,:))+abs(SWIAc.E(jj+1,:)-SWIAc.E(jj,:)))/2;
    end
    SWIAc.dE(1,:) = 1.16*SWIAc.dE(2,:);
    SWIAc.dE(end,:) = 0.87*SWIAc.dE(end-1,:);
    
    SWIAc.E4 = permute(repmat(SWIAc.E,[1 size(SWIAc.ph,1) length(SWIAc.th) length(SWIAc.time)]),[2 3 1 4]);
    SWIAc.dE4 = permute(repmat(SWIAc.dE,[1 size(SWIAc.ph,1) length(SWIAc.th) length(SWIAc.time)]),[2 3 1 4]);
    
    SWIAc.th4 = permute(repmat(SWIAc.th,[1 size(SWIAc.c,1) size(SWIAc.c,3) size(SWIAc.c,4)]),[2 1 3 4]);
    SWIAc.dth = deg2rad(22.5);
    [~,SWIAc.dph] = gradient(SWIAc.ph);
    SWIAc.ph4 = permute(repmat(SWIAc.ph,[1 1 length(SWIAc.th) length(SWIAc.time)]),[1 3 2 4]);
    SWIAc.dph4 = permute(repmat(SWIAc.dph,[1 1 length(SWIAc.th) length(SWIAc.time)]),[1 3 2 4]);
    SWIAc.Omega4 = SWIAc.dph4.*SWIAc.dth.*cos(SWIAc.ph4);
    
    SWIAc.vHE4 = sqrt(2*SWIAc.E4*eV/mH)*100; %[cm/s]
    SWIAc.vAE4 = sqrt(2*2*SWIAc.E4*eV/mA)*100; %[cm/s]
    
    [SWIAc.vHx,SWIAc.vHy,SWIAc.vHz] = sph2cart(SWIAc.th4,SWIAc.ph4,SWIAc.vHE4);
    [SWIAc.vHxn,SWIAc.vHyn,SWIAc.vHzn] = sph2cart(SWIAc.th4,SWIAc.ph4,1);

    [SWIAc.vAx,SWIAc.vAy,SWIAc.vAz] = sph2cart(SWIAc.th4,SWIAc.ph4,SWIAc.vAE4);
    
    %Define blocked angular bins
    [SWIAc.thgLook,SWIAc.phgLook,~] = cart2sph(-SWIAc.vHx,-SWIAc.vHy,-SWIAc.vHz);
    SWIAc.FOV = ones(size(SWIAc.ph4));
    SWIAc.FOV((SWIAc.phgLook < deg2rad(-5)).*(SWIAc.thgLook < deg2rad(-90)) == 1) = 0;
    SWIAc.FOV = (SWIAc.FOV == 1);
    
    SWIAc.fomni = squeeze(nansum(nansum(SWIAc.j.*SWIAc.Omega4,2),1));
    SWIAc.Fomni = squeeze(nansum(nansum(SWIAc.j.*SWIAc.Omega4.*SWIAc.dE4,2),1));
    SWIAc.N = squeeze(nansum(nansum(nansum(SWIAc.j./SWIAc.vHE4.*SWIAc.Omega4.*SWIAc.dE4,3),2),1));
    
    SWIAc.VxSWIA = 1E-5*(1./SWIAc.N).*squeeze(nansum(nansum(nansum(SWIAc.j.*cos(SWIAc.th4).*cos(SWIAc.ph4).*SWIAc.Omega4.*SWIAc.dE4,3),2),1));
    SWIAc.VySWIA = 1E-5*(1./SWIAc.N).*squeeze(nansum(nansum(nansum(SWIAc.j.*sin(SWIAc.th4).*cos(SWIAc.ph4).*SWIAc.Omega4.*SWIAc.dE4,3),2),1));
    SWIAc.VzSWIA = 1E-5*(1./SWIAc.N).*squeeze(nansum(nansum(nansum(SWIAc.j.*sin(SWIAc.ph4).*SWIAc.Omega4.*SWIAc.dE4,3),2),1));
    
    SWIAc.Vx = single(NaN(size(SWIAc.VxSWIA)));
    SWIAc.Vy = single(NaN(size(SWIAc.VySWIA)));
    SWIAc.Vz = single(NaN(size(SWIAc.VzSWIA)));
    for ii = 1:length(SWIAc.time)
        temp = SWIAc.Mswia2msoAll(:,:,ii)*[SWIAc.VxSWIA(ii);SWIAc.VySWIA(ii);SWIAc.VzSWIA(ii)];
        SWIAc.Vx(ii) = temp(1);
        SWIAc.Vy(ii) = temp(2);
        SWIAc.Vz(ii) = temp(3);
    end
    SWIAc.Vabs = sqrt(SWIAc.Vx.^2 + SWIAc.Vy.^2 + SWIAc.Vz.^2);
    
    SWIAc.Pxx = mH*squeeze(nansum(nansum(nansum(SWIAc.j.*10^4.*SWIAc.vHE4/100.*cos(SWIAc.th4).^2.*cos(SWIAc.ph4).^2.*SWIAc.Omega4.*SWIAc.dE4,3),2),1)) - mH*SWIAc.N.*10^6.*(SWIAc.VxSWIA*1000).^2;
    SWIAc.Pyy = mH*squeeze(nansum(nansum(nansum(SWIAc.j.*10^4.*SWIAc.vHE4/100.*sin(SWIAc.th4).^2.*cos(SWIAc.ph4).^2.*SWIAc.Omega4.*SWIAc.dE4,3),2),1)) - mH*SWIAc.N.*10^6.*(SWIAc.VySWIA*1000).^2;
    SWIAc.Pzz = mH*squeeze(nansum(nansum(nansum(SWIAc.j.*10^4.*SWIAc.vHE4/100.*sin(SWIAc.ph4).^2.*SWIAc.Omega4.*SWIAc.dE4,3),2),1)) - mH*SWIAc.N.*10^6.*(SWIAc.VzSWIA*1000).^2;
    
    SWIAc.T = abs(SWIAc.Pxx + SWIAc.Pyy + SWIAc.Pzz)./(3*SWIAc.N*10^6*k); %[K]
    
    SWIAc.Vth = sqrt(2*k*SWIAc.T/mH)/1000; %Thermal speed [km/s]
    
    %Angle relative the SW
    vswn = rotz(-4)*[-1;0;0];
    
    SWIAc.thvsw = acosd(SWIAc.Vx./SWIAc.Vabs.*vswn(1) + SWIAc.Vy./SWIAc.Vabs.*vswn(2) + SWIAc.Vz./SWIAc.Vabs.*vswn(3));
    
    %Solar wind window
    SWIAc.VxnSWIA = SWIAc.VxSWIA./SWIAc.Vabs;
    SWIAc.VynSWIA = SWIAc.VySWIA./SWIAc.Vabs;
    SWIAc.VznSWIA = SWIAc.VzSWIA./SWIAc.Vabs;
    
    SWIAc.isWindow = zeros(size(SWIAc.vHxn)) > 0;
    SWIAc.jE = SWIAc.j.*SWIAc.E4;
    Epeak = NaN(size(SWIAc.time));
    for ii = 1:length(SWIAc.time)
        jE = SWIAc.jE(:,:,:,ii);
        E = SWIAc.E4(:,:,:,ii);
        indjEmax = find(jE == max(jE(:)));
        Epeak(ii) = E(indjEmax(1));
        SWIAc.isWindow(:,:,:,ii) = (acosd(SWIAc.vHxn(:,:,:,ii).*SWIAc.VxnSWIA(ii) + SWIAc.vHyn(:,:,:,ii).*SWIAc.VynSWIA(ii) + SWIAc.vHzn(:,:,:,ii).*SWIAc.VznSWIA(ii)) < 30).*((SWIAc.E4(:,:,:,ii) > Epeak(ii)/3).*(SWIAc.E4(:,:,:,ii) < 3.5*Epeak(ii)) == 1) > 0;
    end
    
    SWIAc.Nwindow = squeeze(nansum(nansum(nansum(SWIAc.j.*SWIAc.isWindow./SWIAc.vHE4.*SWIAc.Omega4.*SWIAc.dE4,3),2),1));
    
    SWIAc.Bx = single(NaN(size(SWIAc.VxSWIA))); %MSO
    SWIAc.By = single(NaN(size(SWIAc.VxSWIA)));
    SWIAc.Bz = single(NaN(size(SWIAc.VxSWIA)));

    SWIAc.BxSWIA = single(NaN(size(SWIAc.VxSWIA))); %SWIA
    SWIAc.BySWIA = single(NaN(size(SWIAc.VxSWIA)));
    SWIAc.BzSWIA = single(NaN(size(SWIAc.VxSWIA)));
    for ii = 1:length(SWIAc.time)
        temp = SWIAc.Mswia2msoAll(:,:,ii)*[SWIAc.VxSWIA(ii);SWIAc.VySWIA(ii);SWIAc.VzSWIA(ii)];
        SWIAc.Vx(ii) = temp(1);
        SWIAc.Vy(ii) = temp(2);
        SWIAc.Vz(ii) = temp(3);
        
        indMAG = (MAG.time > (SWIAc.time(ii) - 2/(24*3600))).*(MAG.time < (SWIAc.time(ii) + 2/(24*3600))) == 1;
        SWIAc.Bx(ii) = mean(MAG.Bx(indMAG));
        SWIAc.By(ii) = mean(MAG.By(indMAG));
        SWIAc.Bz(ii) = mean(MAG.Bz(indMAG));

        tempB = SWIAc.Mmso2swiaAll(:,:,ii)*[SWIAc.Bx(ii);SWIAc.By(ii);SWIAc.Bz(ii)];
        SWIAc.BxSWIA(ii) = tempB(1);
        SWIAc.BySWIA(ii) = tempB(2);
        SWIAc.BzSWIA(ii) = tempB(3);
    end
    SWIAc.Babs = sqrt(SWIAc.Bx.^2 + SWIAc.By.^2 + SWIAc.Bz.^2);
    SWIAc.Bxn = SWIAc.Bx./SWIAc.Babs;
    SWIAc.Byn = SWIAc.By./SWIAc.Babs;
    SWIAc.Bzn = SWIAc.Bz./SWIAc.Babs;
    
    SWIAc.BxnSWIA = SWIAc.BxSWIA./SWIAc.Babs;
    SWIAc.BynSWIA = SWIAc.BySWIA./SWIAc.Babs;
    SWIAc.BznSWIA = SWIAc.BzSWIA./SWIAc.Babs;

    %Find and flag distributions with SW out of the SWIA fine FOV
    [SWIAc.vxn4,SWIAc.vyn4,SWIAc.vzn4] = sph2cart(SWIAc.th4,SWIAc.ph4,1);
    SWIAc.thvsw4 = zeros(size(SWIAc.vxn4));
    SWIAc.flagSWoutofFOV = zeros(size(SWIAc.time)) > 0;
    SWIAc.maxthvsw = zeros(size(SWIAc.time));
    for ii = 1:length(SWIAc.time)
        Mmso2swia = SWIAc.Mmso2swiaAll(:,:,ii);
        vswnSWIA = Mmso2swia*vswn;
        [thvswSWIA,phvswSWIA,~] = cart2sph(vswnSWIA(1),vswnSWIA(2),vswnSWIA(3));
        temp = acosd(SWIAc.vxn4(:,:,:,ii).*vswnSWIA(1) + SWIAc.vyn4(:,:,:,ii).*vswnSWIA(2) + SWIAc.vzn4(:,:,:,ii).*vswnSWIA(3));
        SWIAc.maxthvsw(ii) = min(temp(:));
        
        SWIAc.flagSWoutofFOV(ii) = abs(phvswSWIA) > deg2rad(45);
    end

    %Assign spacecraft potentials to the SWIAc measurements
    if SCpotFileFound == 1 
        [ind,dt] = knnsearch(SCpot.time,SWIAc.time);
        SWIAc.SCpot(dt < 2/(24*3600)) = SCpot.pot(ind(dt < 2/(24*3600)));
        SWIAc.SCpot(dt > 2/(24*3600)) = +8;
    else
        SWIAc.SCpot = 8*ones(size(SWIAc.time));
    end

    %% SWIA Fine
    
    %Remove invalid dists
    indRem = squeeze(max(max(max(SWIAf.c)))) > 1E+5; %Remove indices containing impossibly high counts
    SWIAf.time(indRem) = [];
    SWIAf.E(:,indRem) = [];
    SWIAf.Mmso2swiaAll(:,:,indRem) = [];
    SWIAf.Mswia2msoAll(:,:,indRem) = [];
    SWIAf.atten_state(indRem) = [];
    SWIAf.c(:,:,:,indRem) = [];
    SWIAf.j(:,:,:,indRem) = [];
    SWIAf.ph(:,:,indRem) = [];
    SWIAf.C(:,:,:,indRem) = [];
    SWIAf.xSC(indRem) = [];
    SWIAf.ySC(indRem) = [];
    SWIAf.zSC(indRem) = [];

    SWIAf.dE = single(NaN(size(SWIAf.E)));
    for jj = 2:size(SWIAf.E,1)-1
        SWIAf.dE(jj,:) = (abs(SWIAf.E(jj-1,:)-SWIAf.E(jj,:))+abs(SWIAf.E(jj+1,:)-SWIAf.E(jj,:)))/2;
    end
    SWIAf.dE(1,:) = 1.1*SWIAf.dE(2,:);
    SWIAf.dE(end,:) = 0.9*SWIAf.dE(end-1,:);
    
    SWIAf.E4 = permute(repmat(SWIAf.E,[1 1 size(SWIAf.ph,1) length(SWIAf.th)]),[3 4 1 2]);
    SWIAf.dE4 = permute(repmat(SWIAf.dE,[1 1 size(SWIAf.ph,1) length(SWIAf.th)]),[3 4 1 2]);
    
    SWIAf.th4 = permute(repmat(SWIAf.th,[1 size(SWIAf.c,1) size(SWIAf.c,3) size(SWIAf.c,4)]),[2 1 3 4]);
    SWIAf.dth = deg2rad(4.5);
    [~,SWIAf.dph] = gradient(SWIAf.ph);
    SWIAf.ph4 = permute(repmat(SWIAf.ph,[1 1 1 length(SWIAf.th)]),[1 4 2 3]);
    SWIAf.dph4 = permute(repmat(SWIAf.dph,[1 1 1 length(SWIAf.th)]),[1 4 2 3]);
    SWIAf.Omega4 = SWIAf.dph4.*SWIAf.dth.*cos(SWIAf.ph4);
    
    SWIAf.vHE4 = sqrt(2*SWIAf.E4*eV/mH)*100; %[cm/s]
    SWIAf.vAE4 = sqrt(2*2*SWIAf.E4*eV/mA)*100; %[cm/s]
    
    [SWIAf.vHx,SWIAf.vHy,SWIAf.vHz] = sph2cart(SWIAf.th4,SWIAf.ph4,SWIAf.vHE4);
    [SWIAf.vAx,SWIAf.vAy,SWIAf.vAz] = sph2cart(SWIAf.th4,SWIAf.ph4,SWIAf.vAE4);
    [SWIAf.vxn4,SWIAf.vyn4,SWIAf.vzn4] = sph2cart(SWIAf.th4,SWIAf.ph4,1);

    %Resolution in velocity-space
    dth = deg2rad(4.5);
    dph = deg2rad(3.75);
    dEoverE = 0.145;
    SWIAf.sigvHx = sqrt((SWIAf.vHE4.*(-sin(SWIAf.th4)).*(cos(SWIAf.ph4))).^2.*dth.^2 + (SWIAf.vHE4.*(cos(SWIAf.th4)).*(-sin(SWIAf.ph4))).^2.*dph.^2 + double(1./(sqrt(2*mH*SWIAf.E4*eV)).*cos(SWIAf.th4).*cos(SWIAf.ph4)).^2.*double(dEoverE*SWIAf.E4*eV).^2);
    SWIAf.sigvHy = sqrt((SWIAf.vHE4.*(cos(SWIAf.th4)).*(cos(SWIAf.ph4))).^2.*dth.^2 + (SWIAf.vHE4.*(sin(SWIAf.th4)).*(-sin(SWIAf.ph4))).^2.*dph.^2  + double(1./(sqrt(2*mH*SWIAf.E4*eV)).*sin(SWIAf.th4).*cos(SWIAf.ph4)).^2.*double(dEoverE*SWIAf.E4*eV).^2);
    SWIAf.sigvHz = sqrt((SWIAf.vHE4.*(cos(SWIAf.ph4))).^2.*dph.^2  + double((1./(sqrt(2*mH*SWIAf.E4*eV)).*sin(SWIAf.ph4))).^2.*double(dEoverE*SWIAf.E4*eV).^2);
    
    SWIAf.sigj = SWIAf.j./SWIAf.c.*sqrt(SWIAf.c); %Approximated Poisson uncertainty
    SWIAf.sigjE = SWIAf.j.*SWIAf.E4./SWIAf.c.*sqrt(SWIAf.c); 

    SWIAf.fomni = squeeze(nansum(nansum(SWIAf.j.*SWIAf.Omega4,2),1));
    SWIAf.Fomni = squeeze(nansum(nansum(SWIAf.j.*SWIAf.Omega4.*SWIAf.dE4,2),1));
    SWIAf.N = squeeze(nansum(nansum(nansum(SWIAf.j./SWIAf.vHE4.*SWIAf.Omega4.*SWIAf.dE4,3),2),1));
    
    SWIAf.VxSWIA = 1E-5*(1./SWIAf.N).*squeeze(nansum(nansum(nansum(SWIAf.j.*cos(SWIAf.th4).*cos(SWIAf.ph4).*SWIAf.Omega4.*SWIAf.dE4,3),2),1));
    SWIAf.VySWIA = 1E-5*(1./SWIAf.N).*squeeze(nansum(nansum(nansum(SWIAf.j.*sin(SWIAf.th4).*cos(SWIAf.ph4).*SWIAf.Omega4.*SWIAf.dE4,3),2),1));
    SWIAf.VzSWIA = 1E-5*(1./SWIAf.N).*squeeze(nansum(nansum(nansum(SWIAf.j.*sin(SWIAf.ph4).*SWIAf.Omega4.*SWIAf.dE4,3),2),1));
    
    SWIAf.Vx = single(NaN(size(SWIAf.VxSWIA))); %MSO
    SWIAf.Vy = single(NaN(size(SWIAf.VySWIA)));
    SWIAf.Vz = single(NaN(size(SWIAf.VzSWIA)));

    SWIAf.Bx = single(NaN(size(SWIAf.VxSWIA))); %MSO
    SWIAf.By = single(NaN(size(SWIAf.VxSWIA)));
    SWIAf.Bz = single(NaN(size(SWIAf.VxSWIA)));

    SWIAf.BxSWIA = single(NaN(size(SWIAf.VxSWIA))); %SWIA
    SWIAf.BySWIA = single(NaN(size(SWIAf.VxSWIA)));
    SWIAf.BzSWIA = single(NaN(size(SWIAf.VxSWIA)));
    for ii = 1:length(SWIAf.time)
        temp = SWIAf.Mswia2msoAll(:,:,ii)*[SWIAf.VxSWIA(ii);SWIAf.VySWIA(ii);SWIAf.VzSWIA(ii)];
        SWIAf.Vx(ii) = temp(1);
        SWIAf.Vy(ii) = temp(2);
        SWIAf.Vz(ii) = temp(3);
        
        indMAG = (MAG.time > (SWIAf.time(ii) - 2/(24*3600))).*(MAG.time < (SWIAf.time(ii) + 2/(24*3600))) == 1;
        SWIAf.Bx(ii) = mean(MAG.Bx(indMAG));
        SWIAf.By(ii) = mean(MAG.By(indMAG));
        SWIAf.Bz(ii) = mean(MAG.Bz(indMAG));

        tempB = SWIAf.Mmso2swiaAll(:,:,ii)*[SWIAf.Bx(ii);SWIAf.By(ii);SWIAf.Bz(ii)];
        SWIAf.BxSWIA(ii) = tempB(1);
        SWIAf.BySWIA(ii) = tempB(2);
        SWIAf.BzSWIA(ii) = tempB(3);
    end
    SWIAf.Babs = sqrt(SWIAf.Bx.^2 + SWIAf.By.^2 + SWIAf.Bz.^2);
    SWIAf.Bxn = SWIAf.Bx./SWIAf.Babs;
    SWIAf.Byn = SWIAf.By./SWIAf.Babs;
    SWIAf.Bzn = SWIAf.Bz./SWIAf.Babs;
    
    SWIAf.BxnSWIA = SWIAf.BxSWIA./SWIAf.Babs;
    SWIAf.BynSWIA = SWIAf.BySWIA./SWIAf.Babs;
    SWIAf.BznSWIA = SWIAf.BzSWIA./SWIAf.Babs;
    
    SWIAf.Vabs = sqrt(SWIAf.Vx.^2 + SWIAf.Vy.^2 + SWIAf.Vz.^2);
    
    %Angle relative the typical solar wind vector
    SWIAf.thvsw = acosd(SWIAf.Vx./SWIAf.Vabs.*vswn(1) + SWIAf.Vy./SWIAf.Vabs.*vswn(2) + SWIAf.Vz./SWIAf.Vabs.*vswn(3));
    
    %Find and flag distributions with SW out of the SWIA fine FOV
    SWIAf.vswnSWIA = single(NaN(3,length(SWIAf.time)));
    SWIAf.thvsw4 = single(zeros(size(SWIAf.vxn4)));
    SWIAf.flagSWoutofFOV = zeros(size(SWIAf.time)) > 0;
    SWIAf.maxthvsw = zeros(size(SWIAf.time));
    for ii = 1:length(SWIAf.time)
        Mmso2swia = SWIAf.Mmso2swiaAll(:,:,ii);
        vswnSWIA = Mmso2swia*vswn;
        SWIAf.vswnSWIA(:,ii) = vswnSWIA;
        temp = acosd(SWIAf.vxn4(:,:,:,ii).*vswnSWIA(1) + SWIAf.vyn4(:,:,:,ii).*vswnSWIA(2) + SWIAf.vzn4(:,:,:,ii).*vswnSWIA(3));
        SWIAf.maxthvsw(ii) = min(temp(:));
    end
    SWIAf.flagSWoutofFOV = SWIAf.maxthvsw > 3.5;
    
    %Temperature
    SWIAf.Pxx = mH*squeeze(nansum(nansum(nansum(SWIAf.j.*10^4.*SWIAf.vHE4/100.*cos(SWIAf.th4).^2.*cos(SWIAf.ph4).^2.*SWIAf.Omega4.*SWIAf.dE4,3),2),1)) - mH*SWIAf.N.*10^6.*(SWIAf.VxSWIA*1000).^2;
    SWIAf.Pyy = mH*squeeze(nansum(nansum(nansum(SWIAf.j.*10^4.*SWIAf.vHE4/100.*sin(SWIAf.th4).^2.*cos(SWIAf.ph4).^2.*SWIAf.Omega4.*SWIAf.dE4,3),2),1)) - mH*SWIAf.N.*10^6.*(SWIAf.VySWIA*1000).^2;
    SWIAf.Pzz = mH*squeeze(nansum(nansum(nansum(SWIAf.j.*10^4.*SWIAf.vHE4/100.*sin(SWIAf.ph4).^2.*SWIAf.Omega4.*SWIAf.dE4,3),2),1)) - mH*SWIAf.N.*10^6.*(SWIAf.VzSWIA*1000).^2;
    
    SWIAf.T = abs(SWIAf.Pxx + SWIAf.Pyy + SWIAf.Pzz)./(3*SWIAf.N*10^6*k); %[K]
    
    SWIAf.Vth = sqrt(2*k*SWIAf.T/mH)/1000; %Thermal speed [km/s]

    %Assign spacecraft potentials to the SWIAc measurements
    if SCpotFileFound == 1 && SWIAfFileFound == 1
        [ind,dt] = knnsearch(SCpot.time,SWIAf.time);
        SWIAf.SCpot(dt < 2/(24*3600)) = SCpot.pot(ind(dt < 2/(24*3600)));
        SWIAf.SCpot(dt > 2/(24*3600)) = +8;
    else
        SWIAf.SCpot = 8*ones(size(SWIAf.time));
    end

    %Initialize vectors of manually set parameters for fitting routine
    SWIAf.useRespMat = single(NaN(size(SWIAf.time))); %Force response matrix
    SWIAf.useBeamFastSlow = single(NaN(size(SWIAf.time))); %Force Fast or slow beam (Fast = +1; Slow = -1)
    SWIAf.useBeamBnegBpos = single(NaN(size(SWIAf.time))); %Force Fast or slow beam (Fast = +1; Slow = -1)
    SWIAf.useDVA     = single(NaN(length(SWIAf.time),3));
    SWIAf.useTHcore  = single(NaN(length(SWIAf.time),3)); %Force H+ core temperature ratio constraints
    SWIAf.useTATH    = single(NaN(length(SWIAf.time),3)); %Force He2+/H+ temperature ratio constraints
    SWIAf.useAlphaH  = single(NaN(length(SWIAf.time),3));
    SWIAf.useAlphaA  = single(NaN(length(SWIAf.time),3));
    SWIAf.useElim  = single(NaN(length(SWIAf.time),2));
    SWIAf.useAzLim  = single(NaN(length(SWIAf.time),2));
    SWIAf.useElLim  = single(NaN(length(SWIAf.time),2));

    for ii = 1:size(timeManSWuseindr,1)
        SWIAf.useRespMat((SWIAf.time > timeManSWuseindr(ii,1)).*(SWIAf.time < timeManSWuseindr(ii,2)) == 1) = manSWuseindr(ii);
    end
    for ii = 1:size(timeManSWbeamFast,1)
        SWIAf.useBeamFastSlow((SWIAf.time > timeManSWbeamFast(ii,1)).*(SWIAf.time < timeManSWbeamFast(ii,2)) == 1) = +1;
    end
    for ii = 1:size(timeManSWbeamSlow,1)
        SWIAf.useBeamFastSlow((SWIAf.time > timeManSWbeamSlow(ii,1)).*(SWIAf.time < timeManSWbeamSlow(ii,2)) == 1) = -1;
    end
    for ii = 1:size(timeManSWbeamBneg,1)
        SWIAf.useBeamBnegBpos((SWIAf.time > timeManSWbeamBneg(ii,1)).*(SWIAf.time < timeManSWbeamBneg(ii,2)) == 1) = -1;
    end
    for ii = 1:size(timeManSWbeamBpos,1)
        SWIAf.useBeamBnegBpos((SWIAf.time > timeManSWbeamBpos(ii,1)).*(SWIAf.time < timeManSWbeamBpos(ii,2)) == 1) = +1;
    end
    for ii = 1:size(timeManDVA,1)
        ind = (SWIAf.time > timeManDVA(ii,1)).*(SWIAf.time < timeManDVA(ii,2)) == 1;
        SWIAf.useDVA(ind,:) = repmat(manDVA(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManTHcore,1)
        ind = (SWIAf.time > timeManTHcore(ii,1)).*(SWIAf.time < timeManTHcore(ii,2)) == 1;
        SWIAf.useTHcore(ind,:) = repmat(manTHcore(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManTATH,1)
        ind = (SWIAf.time > timeManTATH(ii,1)).*(SWIAf.time < timeManTATH(ii,2)) == 1;
        SWIAf.useTATH(ind,:) = repmat(manTATH(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManAlphaH,1)
        ind = (SWIAf.time > timeManAlphaH(ii,1)).*(SWIAf.time < timeManAlphaH(ii,2)) == 1;
        SWIAf.useAlphaH(ind,:) = repmat(manAlphaH(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManAlphaA,1)
        ind = (SWIAf.time > timeManAlphaA(ii,1)).*(SWIAf.time < timeManAlphaA(ii,2)) == 1;
        SWIAf.useAlphaA(ind,:) = repmat(manAlphaA(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManElim,1)
        ind = (SWIAf.time > timeManElim(ii,1)).*(SWIAf.time < timeManElim(ii,2)) == 1;
        SWIAf.useElim(ind,:) = repmat(manElim(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManAzLim,1)
        ind = (SWIAf.time > timeManAzLim(ii,1)).*(SWIAf.time < timeManAzLim(ii,2)) == 1;
        SWIAf.useAzLim(ind,:) = repmat(manAzLim(ii,:),[sum(ind),1]);
    end
    for ii = 1:size(timeManElLim,1)
        ind = (SWIAf.time > timeManElLim(ii,1)).*(SWIAf.time < timeManElLim(ii,2)) == 1;
        SWIAf.useElLim(ind,:) = repmat(manElLim(ii,:),[sum(ind),1]);
    end
    
    if true
        disp('   %% Manual fit selections for this orbit (number of distributions):');
        disp(['    %% Response matrix: ',num2str(sum(~isnan(SWIAf.useRespMat)))]);
        disp(['    %% Beam direction: ',num2str(sum(~isnan(SWIAf.useBeamBnegBpos)))]);
        disp(['    %% dVA constraints: ',num2str(sum(~isnan(SWIAf.useDVA)))]);
        disp(['    %% THcore constraints: ',num2str(sum(~isnan(SWIAf.useTHcore)))]);
        disp(['    %% TA/TH constraints: ',num2str(sum(~isnan(SWIAf.useTATH)))]);
        disp(['    %% alphaH constraints: ',num2str(sum(~isnan(SWIAf.useAlphaH)))]);
        disp(['    %% alphaA constraints: ',num2str(sum(~isnan(SWIAf.useAlphaA)))]);
        disp(['    %% Elim intervals: ',num2str(sum(~isnan(SWIAf.useElim)))]);
        disp(['    %% AzLim intervals: ',num2str(sum(~isnan(SWIAf.useAzLim)))]);
        disp(['    %% ElLim intervals: ',num2str(sum(~isnan(SWIAf.useElLim)))]);
        
    end
    
    %% SWEA
    %[SWEA.vex,SWEA.vey,SWEA.vez] = sph2cart(SWEA.thg,SWEA.phg,sqrt(2*SWEA.E*eV/me));
    SWEA.veabs = sqrt(2*SWEA.E*eV/me)*100;
    SWEA.vex = SWEA.veabs.*SWEA.xmso;
    SWEA.vey = SWEA.veabs.*SWEA.ymso;
    SWEA.vez = SWEA.veabs.*SWEA.zmso;
    [SWEA.thg,SWEA.phg,~] = cart2sph(SWEA.xmso,SWEA.ymso,SWEA.zmso);
    SWEA.elevg = repmat(reshape(SWEA.elev,[size(SWEA.elev,1) 1 size(SWEA.elev,2)]),[1 length(SWEA.azim) 1 length(SWEA.time)]);
    SWEA.daz = abs(diff(SWEA.azim)); SWEA.daz(end+1) = SWEA.daz(end);
    SWEA.del = abs(diff(SWEA.elev,1)); SWEA.del(end+1,:) = SWEA.del(1,:);
    SWEA.Omega = repmat(shiftdim(SWEA.daz,-1),[size(SWEA.del,1) 1 size(SWEA.del,2) length(SWEA.time)]).*repmat(reshape(SWEA.del,[size(SWEA.del,1) 1 size(SWEA.del,2)]),[1 length(SWEA.azim) 1 length(SWEA.time)]).*cos(SWEA.elevg);
    SWEA.dE = zeros(size(SWEA.E));
    for ii = 2:size(SWEA.E,3)-1
        SWEA.dE(:,:,ii,:) = (abs(SWEA.E(:,:,ii-1,:)-SWEA.E(:,:,ii,:)) + abs(SWEA.E(:,:,ii,:)-SWEA.E(:,:,ii+1,:)))/2;
    end
    SWEA.j(isnan(SWEA.E)) = NaN;
    SWEA.E2 = squeeze(SWEA.E(3,1,:,:));
    
    SWEA.fomni = squeeze(nansum(nansum(SWEA.j.*SWEA.Omega,2),1));
    SWEA.Fomni = squeeze(nansum(nansum(nansum(SWEA.j.*SWEA.Omega.*SWEA.dE,3),2),1));
    
    SWEA.FomniGtr100eV = single(NaN(size(SWEA.Fomni)));
    SWEA.FomniGtr10eV = single(NaN(size(SWEA.Fomni)));
    for ii = 1:length(SWEA.time)
        indE100 = SWEA.E2(:,ii) > 100;
        indE10 = SWEA.E2(:,ii) > 10;
        SWEA.FomniGtr100eV(ii) = squeeze(nansum(nansum(nansum(SWEA.j(:,:,indE100,ii).*SWEA.Omega(:,:,indE100,ii).*SWEA.dE(:,:,indE100,ii),3),2),1));
        SWEA.FomniGtr10eV(ii) = squeeze(nansum(nansum(nansum(SWEA.j(:,:,indE10,ii).*SWEA.Omega(:,:,indE10,ii).*SWEA.dE(:,:,indE10,ii),3),2),1));
    end
    
    %Total Pressure [nPa]
    SWEA.az4 = repmat(shiftdim(SWEA.azim,1),[6 1 64 length(SWEA.time)]);
    SWEA.Pxx = 1E+9*me*1E+4*squeeze(nansum(nansum(nansum(SWEA.j.*SWEA.veabs.*(cos(SWEA.az4)).^2.*(cos(SWEA.elevg)).^2.*SWEA.Omega.*SWEA.dE,3),2),1));
    SWEA.Pyy = 1E+9*me*1E+4*squeeze(nansum(nansum(nansum(SWEA.j.*SWEA.veabs.*(sin(SWEA.az4)).^2.*(cos(SWEA.elevg)).^2.*SWEA.Omega.*SWEA.dE,3),2),1));
    SWEA.Pzz = 1E+9*me*1E+4*squeeze(nansum(nansum(nansum(SWEA.j.*SWEA.veabs.*(sin(SWEA.elevg)).^2.*SWEA.Omega.*SWEA.dE,3),2),1));
    
    SWEA.P = (SWEA.Pxx + SWEA.Pyy + SWEA.Pzz)/3;
    
    SWEA.N = squeeze(nansum(nansum(SWEA.j./SWEA.veabs.*SWEA.Omega.*SWEA.dE,2),1));
    
    %% MAG 32Hz
    MAG.Babs = sqrt(MAG.Bx.^2 + MAG.By.^2 + MAG.Bz.^2);

    %% Fit the full coarse distributions
    SWIAc.NHfit    = single(NaN(length(SWIAc.time),1));
    SWIAc.NAfit    = single(NaN(length(SWIAc.time),1));
    SWIAc.VHabsFit = single(NaN(length(SWIAc.time),1));
    SWIAc.VAabsFit = single(NaN(length(SWIAc.time),1));
    SWIAc.THfit    = single(NaN(length(SWIAc.time),1));
    SWIAc.TAfit    = single(NaN(length(SWIAc.time),1));
    SWIAc.VHthFit  = single(NaN(length(SWIAc.time),1));
    SWIAc.pfit1d   = single(NaN(7,length(SWIAc.time)));
    SWIAc.JqFit    = single(NaN(48,length(SWIAc.time)));
    SWIAc.fitExitFlag = single(NaN(length(SWIAc.time),1));
    for ii = 1:length(SWIAc.time)

        if SWIAc.N(ii) == 0
            continue;
        end

        clearvars p0 lb ub;
        p0(1) = 1.0*SWIAc.N(ii)*1E+6; %NH [m^-3]
        p0(2) = 0.03*SWIAc.N(ii)*1E+6;%NA [m^-3]
        p0(3) = 0.3*p0(1);
        p0(4) = -SWIAc.Vabs(ii);%Vabs
        p0(5) = SWIAc.T(ii);%TH
        p0(6) = 1;%a |VA|/|VH|
        p0(7) = 4;%b TA/TH

        lb(1) = 0;
        lb(2) = 0;
        lb(3) = 0;
        lb(4) = -SWIAc.Vabs(ii)-50; 
        lb(5) = 0.1*SWIAc.T(ii);
        lb(6) = 0.99;
        lb(7) = 0.8;

        ub(1) = inf;
        ub(2) = 0.2*SWIAc.N(ii)*1E+6;
        ub(3) = 0.5*ub(1);
        ub(4) = -SWIAc.Vabs(ii)+50; if ub(4) > 0; ub(4) = 0; end
        ub(5) = SWIAc.T(ii);
        ub(6) = 1.05;
        ub(7) = 10;
        
        %Magnetic field vector
        E4 = double(SWIAc.E4(:,:,:,ii)); 
        Omega4 = double(SWIAc.Omega4(:,:,:,ii));
        jq = double(SWIAc.j(:,:,:,ii));
        c = double(SWIAc.c(:,:,:,ii));
        C = double(SWIAc.C(:,:,:,ii));
        sigjq = double(C.*sqrt(c+1)./E4);

        c1 = squeeze(sum(sum(c,2),1)); E1 = SWIAc.E;
        Jq = squeeze(sum(sum(jq.*Omega4,2),1));
        sigJq = sqrt(squeeze(sum(sum(sigjq.^2.*Omega4.^2,2),1)));

        %Cylindrically symmetric samples in velocity space
        dth = deg2rad(1);

        [thg,vHabs] = ndgrid(deg2rad(0:1:45)+pi,double(sqrt(2*E1*eV/mH)));
        [vHx,vHr,~] = sph2cart(thg,0,vHabs);

        [thg,vAabs] = ndgrid(deg2rad(0:1:45)+pi,double(sqrt(2*2*E1*eV/mA)));
        [vAx,vAr,~] = sph2cart(thg,0,vAabs);

        %Fitting routine
        clearvars p pfit1d;
        
        Jfun = @(p) squeeze(sum(2*pi*(thg-pi)*dth.*(1E-4.*vHabs.^2./mH.*eV.*p(1).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-p(4)*1000).^2 + vHr.^2)./(2*k*p(5)))... 
                                                   + 2*1E-4.*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(7)*p(5)))^(3/2).*exp(-mA.*((vAx-p(6)*p(4)*1000).^2 + vAr.^2)./(2*k*p(7)*p(5)))),1));%    + 2*1E-4*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(8)*p(6))).^(3/2).*exp(-mA.*((vAx-p(7).*p(3)*1000).^2 + (vAy-p(7).*p(4)*1000).^2 + (vAz-p(7).*p(5)*1000).^2)./(2*k*p(8)*p(6)));
        
        costfun = @(p) sum((Jq - Jfun(p)').^2./sigJq.^2); %Chi-square cost function
        
        opts = optimoptions('fmincon','Display','off');
        [pfit1d,chi2,fitExitFlag,~,~,~,H] = fmincon(costfun,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
       
        chi2r= chi2/(length(Jq)-length(pfit1d));
        
        sigpfit1d = sqrt(diag(inv(H/2)));%.*sqrt(length(Jq)-length(pfit));

        if sum(imag(sigpfit1d) > 0)
            warning('Bad Hessian from 1D fit, solution likely at boundary');
            clearvars sigpfit1d;
            sigpfit1d = inf(1,8);
        end
        SWIAc.NHfit(ii) = pfit1d(1)*1E-6;
        SWIAc.NAfit(ii) = pfit1d(2)*1E-6;
        SWIAc.VHabsFit(ii) = abs(pfit1d(4));
        SWIAc.VAabsFit(ii) = abs(pfit1d(6)*pfit1d(4));
        SWIAc.THfit(ii) = pfit1d(5);
        SWIAc.TAfit(ii) = pfit1d(7)*pfit1d(5);
        
        SWIAc.sigNHfit(ii) = sqrt(sigpfit1d(1)^2+sigpfit1d(3)^2)*1E-6;
        SWIAc.sigNAfit(ii) = sigpfit1d(2)*1E-6;
        SWIAc.sigTHfit(ii) = sigpfit1d(5);
        SWIAc.sigTAfit(ii) = sqrt(pfit1d(2)^2*((sigpfit1d(7)/pfit1d(7))^2 + (sigpfit1d(5)/pfit1d(5))^2));
        
        SWIAc.JqFit(:,ii) = Jfun(pfit1d)';
        
        SWIAc.pfit1d(:,ii) = pfit1d';
        SWIAc.VHthH(ii) = sqrt(2*k*SWIAc.THfit(ii)/mH); %Thermal speed 
        SWIAc.fitExitFlag(ii) = fitExitFlag;

    end

    %% Find solar wind distributions

    %Find the RMS of the B-field and velocity vector variability
    sigBx = single(NaN(size(SWIAc.time)));
    sigBy = single(NaN(size(SWIAc.time)));
    sigBz = single(NaN(size(SWIAc.time)));
    
    SWIAc.BxnAvg = single(NaN(size(SWIAc.time)));
    SWIAc.BynAvg = single(NaN(size(SWIAc.time)));
    SWIAc.BznAvg = single(NaN(size(SWIAc.time)));
    SWIAc.BabsAvg = NaN(size(SWIAc.time));

    sigVx = NaN(size(SWIAc.time));
    sigVy = NaN(size(SWIAc.time));
    sigVz = NaN(size(SWIAc.time));
    VabsAvg = NaN(size(SWIAc.time));
    
    dtSWIAc = gradient(SWIAc.time);
    SWIAc.runStdTrel = single(NaN(size(SWIAc.time)));
    for ii = 1:length(SWIAc.time)
        indSelMag = abs(MAG.time - SWIAc.time(ii)) < 16/(24*3600);
        
        sigBx(ii) = std(MAG.Bx(indSelMag));
        sigBy(ii) = std(MAG.By(indSelMag));
        sigBz(ii) = std(MAG.Bz(indSelMag));
        
        SWIAc.BabsAvg(ii) = mean(MAG.Babs(indSelMag));
        
        indSelSWIAc = abs(SWIAc.time - SWIAc.time(ii)) < 32/(24*3600);
        
        SWIAc.runStdTrel(ii) = nanstd(SWIAc.THfit(indSelSWIAc))./nanmean(SWIAc.THfit(indSelSWIAc));
        
        sigVx(ii) = std(SWIAc.Vx(indSelSWIAc));
        sigVy(ii) = std(SWIAc.Vy(indSelSWIAc));
        sigVz(ii) = std(SWIAc.Vz(indSelSWIAc));
        VabsAvg(ii) = mean(SWIAc.Vabs(indSelSWIAc));
        
    end
    
    SWIAc.sigB = sqrt((sigBx.^2 + sigBy.^2 + sigBz.^2)/3); %RMS deviation in B
    SWIAc.sigV = sqrt((sigVx.^2 + sigVy.^2 + sigVz.^2)/3); %RMS deviation in V
    
    SWIAc.flagSigB = SWIAc.sigB./SWIAc.BabsAvg > 0.15;
    SWIAc.flagSigV = SWIAc.sigV./VabsAvg > 0.05;
    
    %Find ionospheric & strong plume signatures
    if STATICfileFound == 1
        [indSWIAc,dt] = knnsearch(STA.time,SWIAc.time);
        SWIAc.flagHVY = (STA.NO(indSWIAc) + STA.NO2(indSWIAc)) > 0.1;
    else
        SWIAc.flagHVY = zeros(size(SWIAc.time));
    end
    
    %Flag high-energy electrons
    if SWEAfileFound == 1
        [indSWIAc,dt] = knnsearch(SWEA.time,SWIAc.time);
        SWIAc.flagEelec = SWEA.FomniGtr100eV(indSWIAc)./SWEA.Fomni(indSWIAc) > 0.2;
    else
        SWIAc.flagEelec = NaN(size(SWIAc.time));
    end

    %Flag electron pressure >> ion dynamic pressure
    SWIAc.pdyn = 1E+9*0.5*1.1.*mH*1E+6*SWIAc.N.*(SWIAc.Vabs*1000).^2;
    SWIAc.Pe = single(NaN(size(SWIAc.time)));
    if SWEAfileFound == 1
        for ii = 1:length(SWIAc.time)
            indSelSWEA = abs(SWEA.time - SWIAc.time(ii)) < 4/(24*3600);
            SWIAc.Pe(ii)    = mean(SWEA.P(indSelSWEA));
        end
        SWIAc.flagPePdynHigh = SWIAc.Pe./SWIAc.pdyn > 5;
    else
        SWIAc.flagPePdynHigh = false(size(SWIAc.time));
    end
    
    %Flag too slow SW
    SWIAc.flagVslow = SWIAc.Vabs < 230;
    
    %Flag high relative magnetic field variability
    SWIAc.flagVar = SWIAc.runStdTrel > 0.15;
    
    %Flag too hot distributions
    SWIAc.flagHot = zeros(size(SWIAc.time));%SWIAc.VthFine./SWIAc.VabsFine > 0.2;
    
    %Fraction of dist within SW window too low
    SWIAc.flagOutWindw = SWIAc.Nwindow./SWIAc.N < 0.85;
    
    %Deflected SW
    SWIAc.flagDefSW = SWIAc.thvsw > 16;
    
    %Location relative model BS boundaries
    %Bow shock model to identify SW distributions
    [RBSinner,RIMBinner,~] = getBoundaries(SWIAc.xSC/Rm,4,500);
    [~,RIMBnom,~] = getBoundaries(SWIAc.xSC/Rm,2,400);
    [RBSouter,RIMBouter,~] = getBoundaries(SWIAc.xSC/Rm,0.1,270);
    
    SWIAc.RyzSC = sqrt(SWIAc.ySC.^2 + SWIAc.zSC.^2)/Rm;
        
    SWIAc.insBSinner = (SWIAc.RyzSC < RBSinner) == 1;
    SWIAc.insBSouter = (SWIAc.RyzSC < RBSouter) == 1;
    
    SWIAc.insIMBinner = (SWIAc.RyzSC < RIMBinner) == 1;
    SWIAc.insIMBnom = (SWIAc.RyzSC < RIMBnom) == 1;
    SWIAc.insIMBouter = (SWIAc.RyzSC < RIMBouter) == 1;
    
    %Find which incides are inside manually identified solar wind inclusion intervals
    SWIAc.indIncl = zeros(size(SWIAc.time));
    SWIAf.indIncl = zeros(size(SWIAf.time));
    for jj = 1:size(timeManInclSW,1)
        SWIAc.indIncl = SWIAc.indIncl + (SWIAc.time > timeManInclSW(jj,1)).*(SWIAc.time < timeManInclSW(jj,2));
        SWIAf.indIncl = SWIAf.indIncl + (SWIAf.time > timeManInclSW(jj,1)).*(SWIAf.time < timeManInclSW(jj,2));
    end
    SWIAc.indIncl = SWIAc.indIncl > 0;
    SWIAf.indIncl = SWIAf.indIncl > 0;

    %Find which incides are inside manually identified solar wind exclusion intervals
    SWIAc.indExcl = zeros(size(SWIAc.time));
    SWIAf.indExcl = zeros(size(SWIAf.time));
    for jj = 1:size(timeManExclSW,1)
        SWIAc.indExcl = SWIAc.indExcl + (SWIAc.time > timeManExclSW(jj,1)).*(SWIAc.time < timeManExclSW(jj,2));
        SWIAf.indExcl = SWIAf.indExcl + (SWIAf.time > timeManExclSW(jj,1)).*(SWIAf.time < timeManExclSW(jj,2));
    end
    SWIAc.indExcl = SWIAc.indExcl > 0;
    SWIAf.indExcl = SWIAf.indExcl > 0;
    
    %Mark plasma domains
    SWIAc.isIono = SWIAc.flagHVY;
    
    SWIAc.isSW = (SWIAc.flagSigB + SWIAc.flagSigV + SWIAc.flagEelec + SWIAc.flagHVY + SWIAc.flagVslow + SWIAc.flagDefSW + SWIAc.flagPePdynHigh + SWIAc.flagVar + SWIAc.flagHot + SWIAc.flagOutWindw) == 0;
    SWIAc.isSW(SWIAc.insBSinner) = false;
    SWIAc.isSW(SWIAc.insBSouter == 0) = true;
    SWIAc.isSW((SWIAc.sigB./SWIAc.BabsAvg > 0.3).*SWIAc.flagVar == 1) = false; % Exclude particularly strong foreshock signatures
    SWIAc.isSW((SWIAc.sigV./SWIAc.Vabs > 0.05) == 1) = false; % Exclude particularly strong foreshock signatures
    SWIAc.isSW(SWIAc.N > 200) = false;
    SWIAc.isSW(SWIAc.indIncl) = true;
    SWIAc.isSW(SWIAc.flagSWoutofFOV) = false;
    SWIAc.isSW(SWIAc.indExcl) = false;

    SWIAc.isIMF = (SWIAc.flagSigB + SWIAc.flagEelec + SWIAc.flagHVY + SWIAc.flagVslow + SWIAc.flagDefSW + SWIAc.flagPePdynHigh + SWIAc.flagVar + SWIAc.flagHot) == 0;
    SWIAc.isIMF(SWIAc.insBSinner) = false;
    SWIAc.isIMF(SWIAc.insBSouter == 0) = true;
    SWIAc.isIMF((SWIAc.sigB./SWIAc.BabsAvg > 0.3).*SWIAc.flagVar == 1) = false;%Exclude particularly strong foreshock signatures
    SWIAc.isIMF(SWIAc.N > 200) = false;
    SWIAc.isIMF(SWIAc.indIncl) = true;
    SWIAc.isIMF(SWIAc.indExcl) = false;
    
    SWIAc.isFS = SWIAc.flagSigB.*SWIAc.isIMF == 1; %Foreshock
    
    SWIAc.isDraped = (SWIAc.flagSigB + SWIAc.flagEelec + SWIAc.flagDefSW) > 0;
    SWIAc.isDraped(SWIAc.insIMBnom) = false;
    SWIAc.isDraped(SWIAc.isSW) = false;
    SWIAc.isDraped(SWIAc.isIono == 1) = false;
    SWIAc.isDraped(SWIAc.indIncl) = false;

    %% Create a joint time-series based on Coarse and Fine distributions
    clearvars SWIA; 
    numelsE = 200;
    SWIAc.dtime = NaN(size(SWIAc.time)); %Time difference to nearest neighbor distribution
    for ii = 2:length(SWIAc.time)-1
        SWIAc.dtime(ii) = min([SWIAc.time(ii+1) - SWIAc.time(ii), SWIAc.time(ii) - SWIAc.time(ii-1)]);
    end
    SWIAc.dtime(1) = SWIAc.time(2)-SWIAc.time(1);
    SWIAc.dtime(end) = SWIAc.time(end)-SWIAc.time(end-1);
    SWIAc.dtime = abs(SWIAc.dtime);
    
    SWIA.time      = double(zeros(0,1)); % Datenum time
    SWIA.NH        = single(zeros(0,1)); % H+ bulk density [cm^-3]
    SWIA.NHcore    = single(zeros(0,1)); % H+ core density [cm^-3]
    SWIA.NHbeam    = single(zeros(0,1)); % H+ beam density [cm^-3]
    SWIA.NA        = single(zeros(0,1)); % He2+ (total) density [cm^-3]

    SWIA.VHx       = single(zeros(0,1)); % Bulk H+ velocity x-component in MSO reference frame [km/s]
    SWIA.VHy       = single(zeros(0,1)); % y-component
    SWIA.VHz       = single(zeros(0,1)); % z-component
    SWIA.VHabs     = single(zeros(0,1)); % H+ speed

    SWIA.VHcorex   = single(zeros(0,1)); % Core H+ velocity x-component in MSO reference frame [km/s]
    SWIA.VHcorey   = single(zeros(0,1)); % y-component
    SWIA.VHcorez   = single(zeros(0,1)); % z-component

    SWIA.dVHbeamx  = single(zeros(0,1)); % Core H+ differential velocity x-component in MSO reference frame [km/s]
    SWIA.dVHbeamy  = single(zeros(0,1)); % y-component
    SWIA.dVHbeamz  = single(zeros(0,1)); % z-component
    
    SWIA.VAx       = single(zeros(0,1)); % He++ velocity x-component in MSO reference frame [km/s]
    SWIA.VAy       = single(zeros(0,1));
    SWIA.VAz       = single(zeros(0,1));
    SWIA.VAabs     = single(zeros(0,1)); % He++ speed

    SWIA.THcore    = single(zeros(0,1)); % H+ core temperature [K]
    SWIA.TA        = single(zeros(0,1)); % He++ Temperature [K]
    SWIA.TH        = single(zeros(0,1)); % H+ bulk temperature (treating all H as single distribution) [K]
    SWIA.THperp    = single(zeros(0,1)); % H+ bulk perpendicular temperature (treating all H as single distribution) [K]
    SWIA.THpara    = single(zeros(0,1)); % H+ bulk paralell temperature (treating all H as single distribution) [K]
    SWIA.alphaHcore = single(zeros(0,1)); % H+ core temperature anisotropy factor
    SWIA.alphaA    = single(zeros(0,1)); % Temperature anisotropy factor
    SWIA.kappa     = single(zeros(0,1)); % Kappa (thermilization) factor
    
    SWIA.Bx        = single(zeros(0,1));
    SWIA.By        = single(zeros(0,1));
    SWIA.Bz        = single(zeros(0,1));

    SWIA.BxSWIA    = single(zeros(0,1));
    SWIA.BySWIA    = single(zeros(0,1));
    SWIA.BzSWIA    = single(zeros(0,1));

    SWIA.vnxSWIA   = single(zeros(0,1)); % H+ beam direction differential motion in SWIA reference frame
    SWIA.vnySWIA   = single(zeros(0,1)); % H+ beam direction differential motion in SWIA reference frame
    SWIA.vnzSWIA   = single(zeros(0,1)); % H+ beam direction differential motion in SWIA reference frame
    
    SWIA.Valf      = single(zeros(0,1)); %Alfv?n wave speed

    SWIA.En          = single(NaN(0,numelsE)); %Natural (deconvolved) energy spectra 
    SWIA.jEnCoreOmni = single(NaN(0,numelsE));
    SWIA.jEnBeamOmni = single(NaN(0,numelsE));
    SWIA.jEnAomni    = single(NaN(0,numelsE));
    SWIA.jEnRingOmni = single(NaN(0,numelsE));

    SWIA.sigNH     = single(zeros(0,1));
    SWIA.sigNHcore = single(zeros(0,1));
    SWIA.sigNHbeam = single(zeros(0,1));
    SWIA.sigNA     = single(zeros(0,1));

    SWIA.sigVHx    = single(zeros(0,1));
    SWIA.sigVHy    = single(zeros(0,1));
    SWIA.sigVHz    = single(zeros(0,1));

    SWIA.sigVHcorex = single(zeros(0,1));
    SWIA.sigVHcorey = single(zeros(0,1));
    SWIA.sigVHcorez = single(zeros(0,1));

    SWIA.sigdVHbeamx = single(zeros(0,1));
    SWIA.sigdVHbeamy = single(zeros(0,1));
    SWIA.sigdVHbeamz = single(zeros(0,1));

    SWIA.sigVAz    = single(zeros(0,1));
    SWIA.sigVAx    = single(zeros(0,1));
    SWIA.sigVAy    = single(zeros(0,1));

    SWIA.sigTHcore = single(zeros(0,1));
    SWIA.sigTA     = single(zeros(0,1));
    SWIA.sigAlphaHcore = single(zeros(0,1));
    SWIA.sigAlphaA = single(zeros(0,1));
    SWIA.sigKappa  = single(zeros(0,1));

    SWIA.Jq        = single(zeros(0,48));
    SWIA.Eq        = single(zeros(0,48));
    SWIA.pfit      = single(zeros(0,15));
    SWIA.frac95    = single(zeros(0,1));
    SWIA.chi2rApprox = single(zeros(0,1));
    
    SWIA.indrSel    = single(zeros(0,1));
    SWIA.atten_state    = single(zeros(0,1));

    SWIA.isSW      = false(0,1);
    SWIA.isIMF     = false(0,1);
    SWIA.isFS      = false(0,1);
    SWIA.isDraped  = false(0,1);
    SWIA.isIono    = false(0,1);
    SWIA.isFine    = false(0,1);
    SWIA.isSvy     = false(0,1);
    SWIA.fitExitFlag  = single(zeros(0,1));
    SWIA.Mmso2swia = zeros(0,3,3);
    SWIA.Mswia2mso = zeros(0,3,3);
    
    SWIA.xSC  = single(zeros(0,1));
    SWIA.ySC  = single(zeros(0,1));
    SWIA.zSC  = single(zeros(0,1));

    SWIA.dtRunTime = single(zeros(0,1));

    indsFineLog = single(zeros(0,1));
    
    %% Main fitting routine
    indsFineProc = single(zeros(1,0));
    for ii = 1:length(SWIAc.time)
        if SWIAc.isSW(ii)
            indsFine = find(abs(SWIAf.time-SWIAc.time(ii)) < SWIAc.dtime(ii));
            indsFineIncl = true(size(indsFine)); %Log the fine indices that have been included
            
            if ~isempty(indsFine)
                jj = 0;
                while true %Fit the fine distributions
                    tic;
                    jj = jj+1;
                    
                    if jj > length(indsFine)
                        break;
                    end

                    indf = indsFine(jj);

                    %Check if fine distribution is likely missing the bulk distribution
                    if (SWIAf.grouping(indf) == 0) && (SWIAf.N(indf) < 0.8*SWIAc.Nwindow(ii))
                        indsFineIncl(jj) = false;
                        continue;
                    end
                    
                    %Check if this fine distribution has already been processed
                    if ismember(indf,indsFineProc)
                        continue;
                    end
                    indsFineProc = cat(2,indsFineProc,indf);


                    %Initial 1D estimate
                    clearvars p0 lb ub;
                    p0(1) = 1.0*SWIAf.N(indf)*1E+6; %NH [m^-3]
                    p0(2) = 0.03*SWIAf.N(indf)*1E+6;%NA [m^-3]
                    p0(3) = 0.3*p0(1);
                    p0(4) = -SWIAf.Vabs(indf);%Vabs
                    p0(5) = SWIAf.T(indf);%TH
                    p0(6) = 1;%a |VA|/|VH|
                    p0(7) = 4;%b TA/TH
            
                    lb(1) = 0.5*SWIAf.N(indf)*1E+6;
                    lb(2) = 0;
                    lb(3) = 0;
                    lb(4) = -SWIAf.Vabs(indf)-50;
                    lb(5) = 0.5*SWIAf.T(indf);
                    lb(6) = 0.99;
                    lb(7) = 0.8;
            
                    ub(1) = 1.5*SWIAf.N(indf)*1E+6;
                    ub(2) = 0.2*SWIAf.N(indf)*1E+6;
                    ub(3) = 0.5*ub(1);
                    ub(4) = -SWIAf.Vabs(indf)+50;
                    ub(5) = 10*SWIAf.T(indf);
                    ub(6) = 1.05;
                    ub(7) = 10;

                    %Skip empty distributions
                    if sum(isnan([p0,lb,ub])) > 0
                        continue;
                    end
                    
                    %If magnetic field measurement is not available for
                    %this SWIA Fine measurement, use coarse
                    if ~isnan(SWIAf.Babs(indf))

                        %Magnetic field vector
                        Bxn = double(SWIAf.BxnSWIA(indf)); 
                        Byn = double(SWIAf.BynSWIA(indf)); 
                        Bzn = double(SWIAf.BznSWIA(indf));
                        Babs = double(SWIAf.Babs(indf));
                        
                        foundBvector = true;

                    elseif ~isnan(SWIAc.Babs(ii))

                        SWIAf.BxnSWIA(indf) = SWIAc.BxnSWIA(ii);
                        SWIAf.BynSWIA(indf) = SWIAc.BynSWIA(ii);
                        SWIAf.BznSWIA(indf) = SWIAc.BznSWIA(ii);
                        SWIAf.Babs(indf) = SWIAc.Babs(ii);

                        %Magnetic field vector
                        Bxn = double(SWIAf.BxnSWIA(indf)); 
                        Byn = double(SWIAf.BynSWIA(indf)); 
                        Bzn = double(SWIAf.BznSWIA(indf));
                        Babs = double(SWIAf.Babs(indf));
                        
                        foundBvector = true;

                    else %If none available at all, skip it
                        continue;
                    end

                    %Extract selected distribution
                    E4 = double(SWIAf.E4(:,:,:,indf)); 
                    Omega4 = double(SWIAf.Omega4(:,:,:,indf));
                    jq = double(SWIAf.j(:,:,:,indf));
                    c = double(SWIAf.c(:,:,:,indf));
                    C = double(SWIAf.C(:,:,:,indf));
                    sigjq = double(C.*sqrt(c+1)./E4);
                    
                    %1D parameters
                    c1 = squeeze(sum(sum(c,2),1)); E1 = SWIAf.E(:,indf);
                    Jq = squeeze(sum(sum(jq.*Omega4,2),1));
                    sigJq = sqrt(squeeze(sum(sum(sigjq.^2.*Omega4.^2,2),1)));

                    %Cylindrically symmetric samples in velocity space
                    dth = deg2rad(1);
            
                    [thg,vHabs] = ndgrid(deg2rad(0.5:1:45.5)+pi,double(sqrt(2*E1*eV/mH)));
                    [vHx,vHr,~] = sph2cart(thg,0,vHabs);
            
                    [thg,vAabs] = ndgrid(deg2rad(0.5:1:45.5)+pi,double(sqrt(2*2*E1*eV/mA)));
                    [vAx,vAr,~] = sph2cart(thg,0,vAabs);
            
                    [~,E2] = ndgrid(deg2rad(0.5:1:45.5)+pi,double(E1));

                    %Fitting routine
                    clearvars p pfit1d;
                    
                    JfunSlowBeam = @(p) squeeze(sum(2*pi*(thg-pi)*dth.*(1E-4.*vHabs.^2./mH.*eV.*p(1).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-p(4)*1000).^2 + vHr.^2)./(2*k*p(5)))... % Core
                                                                      + 1E-4.*vHabs.^2./mH.*eV.*p(3).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-(p(4)*1000+abs(Bxn)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2))))).^2 + vHr.^2)./(2*k*p(5)))... % Beam
                                                                    + 2*1E-4.*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(7)*p(5)))^(3/2).*exp(-mA.*((vAx-p(6)*p(4)*1000).^2 + vAr.^2)./(2*k*p(7)*p(5)))),1));% Alphas
                            
                    JfunFastBeam = @(p) squeeze(sum(2*pi*(thg-pi)*dth.*(1E-4.*vHabs.^2./mH.*eV.*p(1).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-p(4)*1000).^2 + vHr.^2)./(2*k*p(5)))... 
                                                                      + 1E-4.*vHabs.^2./mH.*eV.*p(3).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-(p(4)*1000-abs(Bxn)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2))))).^2 + vHr.^2)./(2*k*p(5)))...
                                                                    + 2*1E-4.*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(7)*p(5)))^(3/2).*exp(-mA.*((vAx-p(6)*p(4)*1000).^2 + vAr.^2)./(2*k*p(7)*p(5)))),1));
                    
                    %Cost-functions for fast and slow beam models
                    costfunSlowBeam  = @(p) sum((Jq - JfunSlowBeam(p)').^2./sigJq.^2); %Chi-square cost function
                    costfunFastBeam  = @(p) sum((Jq - JfunFastBeam(p)').^2./sigJq.^2); %Chi-square cost function
                    
                    %Fitting
                    opts = optimoptions('fmincon','Display','off');
                    [pfit1dSlow,chi2Slow,fitExitFlagSlow,~,~,~,Hslow] = fmincon(costfunSlowBeam,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
                    [pfit1dFast,chi2Fast,fitExitFlagFast,~,~,~,Hfast] = fmincon(costfunFastBeam,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
                    
                    %Reduced chi-square values
                    chi2rSlow = chi2Slow/(length(Jq)-length(pfit1dSlow));
                    chi2rFast = chi2Fast/(length(Jq)-length(pfit1dFast));

                    %Select best-fitting or manually selected model
                    if SWIAf.useBeamFastSlow(indf) == +1 %Fast beam
                        pfit1d = pfit1dFast;
                        chi2r = chi2rFast;
                        fitExitFlag = fitExitFlagFast;
                        H = Hfast;
                        Jfun = JfunFastBeam;
                        isFastBeam = 1;
                    elseif SWIAf.useBeamFastSlow(indf) == -1 %Slow beam
                        pfit1d = pfit1dSlow;
                        chi2r = chi2rSlow;
                        fitExitFlag = fitExitFlagSlow;
                        H = Hslow;
                        Jfun = JfunSlowBeam;
                        isFastBeam = 0;
                    else %No manual selection
                        if chi2rFast < chi2rSlow
                            pfit1d = pfit1dFast;
                            chi2r = chi2rFast;
                            fitExitFlag = fitExitFlagFast;
                            H = Hfast;
                            Jfun = JfunFastBeam;
                            isFastBeam = 1;
                        else
                            pfit1d = pfit1dSlow;
                            chi2r = chi2rSlow;
                            fitExitFlag = fitExitFlagSlow;
                            H = Hslow;
                            Jfun = JfunSlowBeam;
                            isFastBeam = 0;
                        end
                    end
                    
                    %Derive uncertainties from the diagonal elements of the Hessian
                    sigpfit1d = sqrt(diag(inv(H/2)));
        
                    if sum(imag(sigpfit1d) > 0)
                        warning('Bad Hessian from 1D fit, solution likely at boundary');
                        clearvars sigpfit1d;
                        sigpfit1d(1) = pfit1d(1);
                        sigpfit1d(2) = pfit1d(2);
                        sigpfit1d(3) = pfit1d(3);
                        sigpfit1d(4) = 100;
                        sigpfit1d(5) = 1E+5;
                        sigpfit1d(6) = 0.1;
                        sigpfit1d(7) = pfit1d(7);
                    end
                    
                    pfit1d(4) = abs(pfit1d(4)); %Speed is abs(Vx)

                    % Fit full 3D distribution
                    clearvars p0 lb ub tempH tempA;
                    p0(1) = pfit1d(1); %H+ core density NHcore [m^-3]
                    p0(2) = pfit1d(2); %He++ density NA [m^-3]
                    p0(3) = 0.1*pfit1d(1); %H+ beam density NHbeam [m^-3]
                    p0(4) = pfit1d(4)*(SWIAf.VxSWIA(indf)./SWIAf.Vabs(indf));
                    p0(5) = pfit1d(4)*(SWIAf.VySWIA(indf)./SWIAf.Vabs(indf));
                    p0(6) = pfit1d(4)*(SWIAf.VzSWIA(indf)./SWIAf.Vabs(indf));
                    p0(7) = pfit1d(5)/1.5;%THcore
                    p0(8) = 0;%He++ differential speed along B in Alfven wave speeds
                    p0(9) = 1;%TA/THcore
                    p0(10) = 1; %Alpha H+ (THcore_perp/THcore_para)
                    p0(11) = 1; %Alpha for He++ (TA_perp/TA_para)
                    p0(12) = 4.5; %Kappa factor
                    p0(13) = 1; %Beam differential speed along B in Alfven wave speeds
                    p0(14) = 0.5*1E+6; %Ring total density
                    p0(15) = 0.05; %Average background noise count
            
                    lb(1) = 0;
                    lb(2) = 0;
                    lb(3) = 0;
                    lb(4) = p0(4)-150;
                    lb(5) = p0(5)-150;
                    lb(6) = p0(6)-150;
                    lb(7) = 5E+3;
                    lb(8) = -1.2;
                    lb(9) = pfit1d(7)-3*sigpfit1d(7); if lb(9) < 0.5; lb(9) = 0.5; end
                    lb(10) = 0.1;
                    lb(11) = 0.1;
                    lb(12) = 3/2*1.05;
                    lb(13) = 0.95;
                    lb(14) = 0;
                    lb(15) = 1E-9;
            
                    ub(1) = inf;
                    ub(2) = p0(1);
                    ub(3) = 0.5*p0(1);
                    ub(4) = p0(4)+100;
                    ub(5) = p0(5)+100;
                    ub(6) = p0(6)+100;
                    ub(7) = pfit1d(5)+3*sigpfit1d(5);
                    ub(8) = +1.2; 
                    ub(9) = pfit1d(7)+3*sigpfit1d(7); if ub(9) > 20; ub(9) = 20; end
                    ub(10) = 10;
                    ub(11) = 10;
                    ub(12) = 162;
                    ub(13) = 4;
                    ub(14) = 100*1E+6;
                    ub(15) = 1;

                    %Extract the distribution
                    E4 = double(SWIAf.E4(:,:,:,indf)); E1 = squeeze(E4(1,1,:));
                    dE4 = double(SWIAf.dE4(:,:,:,indf));
                    Omega4 = double(SWIAf.Omega4(:,:,:,indf));
                    th4 = double(SWIAf.th4(:,:,:,indf)); ph4 = double(SWIAf.ph4(:,:,:,indf));
                    vHx = double(SWIAf.vHx(:,:,:,indf)/100);    vAx = double(SWIAf.vAx(:,:,:,indf)/100);
                    vHy = double(SWIAf.vHy(:,:,:,indf)/100);    vAy = double(SWIAf.vAy(:,:,:,indf)/100);
                    vHz = double(SWIAf.vHz(:,:,:,indf)/100);    vAz = double(SWIAf.vAz(:,:,:,indf)/100);
                    vHabs = double(SWIAf.vHE4(:,:,:,indf)/100); vAabs = double(SWIAf.vAE4(:,:,:,indf)/100);
                    
                    jq = double(SWIAf.j(:,:,:,indf));
                    c = double(SWIAf.c(:,:,:,indf)); c1 = squeeze(sum(sum(c,2),1));
                    C = double(SWIAf.C(:,:,:,indf));
                    sigjq = double(C.*sqrt(c+1)./E4);
                    
                    %Skip distribution if practically empty
                    if sum(c(:) > 0) <= 1
                        continue;
                    end
                    
                    %Use specially assigned parameter constraints
                    if ~isnan(SWIAf.useDVA(indf,1))
                        lb(8) = SWIAf.useDVA(indf,1);
                        p0(8) = SWIAf.useDVA(indf,2);
                        ub(8) = SWIAf.useDVA(indf,3);
                    end
                    if ~isnan(SWIAf.useTHcore(indf,1))
                        lb(7) = SWIAf.useTHcore(indf,1);
                        p0(7) = SWIAf.useTHcore(indf,2);
                        ub(7) = SWIAf.useTHcore(indf,3);
                    end
                    if ~isnan(SWIAf.useTATH(indf,1))
                        lb(9) = SWIAf.useTATH(indf,1);
                        p0(9) = SWIAf.useTATH(indf,2);
                        ub(9) = SWIAf.useTATH(indf,3);
                    end
                    if ~isnan(SWIAf.useAlphaH(indf,1))
                        lb(10) = SWIAf.useAlphaH(indf,1);
                        p0(10) = SWIAf.useAlphaH(indf,2);
                        ub(10) = SWIAf.useAlphaH(indf,3);
                    end
                    if ~isnan(SWIAf.useAlphaA(indf,1))
                        lb(11) = SWIAf.useAlphaA(indf,1);
                        p0(11) = SWIAf.useAlphaA(indf,2);
                        ub(11) = SWIAf.useAlphaA(indf,3);
                    end

                    %Remove zeroed elements if survey data
                    if SWIAf.grouping(indf) == 1
                        vHx = vHx(3:10,3:8,9:40);      vAx = vAx(3:10,3:8,9:40);
                        vHy = vHy(3:10,3:8,9:40);      vAy = vAy(3:10,3:8,9:40);
                        vHz = vHz(3:10,3:8,9:40);      vAz = vAz(3:10,3:8,9:40);
                        vHabs = vHabs(3:10,3:8,9:40);  vAabs = vAabs(3:10,3:8,9:40);
            
                        jq = jq(3:10,3:8,9:40);
                        c = c(3:10,3:8,9:40);
                        C = C(3:10,3:8,9:40);
                        sigjq = sigjq(3:10,3:8,9:40);
                        ph4 = ph4(3:10,3:8,9:40);
                        th4 = th4(3:10,3:8,9:40);
                        E4  =  E4(3:10,3:8,9:40);
                        dE4 = dE4(3:10,3:8,9:40);
                        Omega4 = Omega4(3:10,3:8,9:40);
                        E1 = E1(9:40);
                        c1 = c1(9:40);
                    end
                    
                    %Select suitable response function (attenuator dependent)
                    [VthSWIA,VphSWIA,~] = cart2sph(SWIAf.VxSWIA(indf),SWIAf.VySWIA(indf),SWIAf.VzSWIA(indf));

                    Vth1d = sqrt(2*pfit1d(5)*k/mH); %Themal speed from 1D fit
                    thWidth = atand(Vth1d./(1000*pfit1d(4))); %Estimated thermal angular width
                    dEHthEb1d = pfit1d(5)*k/(mH*(1000*pfit1d(4)).^2/2); %Estimated relative thermal energy width (dE/E)
                    dEAthEb1d = pfit1d(7)*pfit1d(5)*k/(mA*(1000*pfit1d(4)).^2/2)/2; %Estimated relative thermal energy width (dE/E)
                    if  abs(rad2deg(VphSWIA)) > 30 && 20*min([dEHthEb1d dEAthEb1d]) < (0.2/7.8)%<0.025
                        indrSel = 4;
                    elseif 20*min([dEHthEb1d dEAthEb1d]) < (0.2/7.8)
                        indrSel = 3;
                    elseif 15*min([dEHthEb1d dEAthEb1d]) < (0.4/7.8) || abs(rad2deg(VphSWIA)) > 25 %<0.05
                        indrSel = 2;
                    elseif 10*min([dEHthEb1d dEAthEb1d]) < (1/7.8) && abs(rad2deg(VphSWIA)) < 25 %<0.13
                        indrSel = 1;
                    else
                        indrSel = 2;
                    end

                    %Use any pre-assigned response matricies
                    if ~isnan(SWIAf.useRespMat(indf))
                        indrSel = SWIAf.useRespMat(indf);
                    end
                    
                    %Assign the response function to be used
                    if indrSel > 0
                        if SWIAf.atten_state(indf) == 1 %Attenuator open
                            SWIAr = SWIAro{indrSel};
                        elseif SWIAf.atten_state(indf) == 2 %Attenuator closed
                            SWIAr = SWIArc{indrSel};
                        end
                    end

                    %Use any potential manually assigned croppings
                    if ~isnan(sum(SWIAf.useElim(indf,:)))
                        indEsel = (E1 > SWIAf.useElim(indf,1)).*(E1 < SWIAf.useElim(indf,2)) == 1;
                    else
                        indEsel = true(size(c,3),1);
                    end
                    if ~isnan(sum(SWIAf.useAzLim(indf,:)))
                        indAzSel = (squeeze(rad2deg(th4(1,:,1))) > SWIAf.useAzLim(indf,1)).*(squeeze(rad2deg(th4(1,:,1))) < SWIAf.useAzLim(indf,2)) == 1;
                    else
                        indAzSel = true(size(c,2),1);
                    end
                    if ~isnan(sum(SWIAf.useElLim(indf,:)))
                        indElSel = (squeeze(rad2deg(ph4(:,1,find(indEsel,1,'last')))) > SWIAf.useElLim(indf,1)).*(squeeze(rad2deg(ph4(:,1,find(indEsel,1,'last')))) < SWIAf.useElLim(indf,2)) == 1;
                    else
                        indElSel = true(size(c,1),1);
                    end
                    cthres = 1.5; %Count threshold
                    if (sum(indEsel) < length(E1)) || (sum(indAzSel) < size(c,2)) %Use manual cropping if defined in either dimension
                        isCropped = true;

                        vHx = vHx(indElSel,indAzSel,indEsel);      vAx = vAx(indElSel,indAzSel,indEsel);
                        vHy = vHy(indElSel,indAzSel,indEsel);      vAy = vAy(indElSel,indAzSel,indEsel);
                        vHz = vHz(indElSel,indAzSel,indEsel);      vAz = vAz(indElSel,indAzSel,indEsel);
                        vHabs = vHabs(indElSel,indAzSel,indEsel);  vAabs = vAabs(indElSel,indAzSel,indEsel);
            
                        jq     =    jq(indElSel,indAzSel,indEsel);
                        c      =     c(indElSel,indAzSel,indEsel);
                        C      =     C(indElSel,indAzSel,indEsel);
                        sigjq  = sigjq(indElSel,indAzSel,indEsel);
                        ph4    =   ph4(indElSel,indAzSel,indEsel);
                        th4    =   th4(indElSel,indAzSel,indEsel);
                        E4     =    E4(indElSel,indAzSel,indEsel);
                        dE4    =   dE4(indElSel,indAzSel,indEsel);
                        Omega4 = Omega4(indElSel,indAzSel,indEsel);
                        E1     = E1(indEsel);
                        c1     = c1(indEsel);

                        indLim1 = [find(indElSel,1,'first') find(indElSel,1,'last')];%[1 size(c,1)];
                        indLim2 = [find(indAzSel,1,'first') find(indAzSel,1,'last')];
                        indLim3 = [find(indEsel,1,'first') find(indEsel,1,'last')];

                    elseif indrSel > 0 && max(c(:)) > cthres %Automatically crop distributions if using high-resolution response function
                    
                        if indrSel == length(SWIArc)
                            nmarg = [1 1 3]; %Margins around high-count area (# indicies; az,el,E/q)
                        elseif size(SWIArc{3}.ksAll,1) > 50
                            nmarg = [1 1 3]; %Margins around high-count area
                        else 
                            nmarg = [1 2 3]; %Margins around high-count area 
                        end
                        isCropped = false;
                        if size(SWIAr.ksAll,1) > 8
                            isCropped = true;
    
                            maxc1 = repmat(max(c,[],1),[size(c,1) 1 1]); 
                            maxc2 = repmat(max(c,[],2),[1 size(c,2) 1]); 
                            maxc3 = repmat(max(c,[],3),[1 1 size(c,3)]);
                            
                            %Cube of max counts over threshold
                            indSel = ((maxc1 > cthres).*(maxc2 > cthres).*(maxc3 > cthres)) > 0;
    
                            [inds1,inds2,inds3] = ind2sub(size(indSel),find(indSel)); 
    
                            indLim1 = [min(inds1) max(inds1)] + nmarg(1)*[-1 +1];
                            indLim2 = [min(inds2) max(inds2)] + nmarg(2)*[-1 +1];
                            indLim3 = [min(inds3) max(inds3)] + nmarg(3)*[-1 +1];
    
                            indSel1 = max([indLim1(1) 1]):min([indLim1(2) size(c,1)]);
                            indSel2 = max([indLim2(1) 1]):min([indLim2(2) size(c,2)]);
                            indSel3 = max([indLim3(1) 1]):min([indLim3(2) size(c,3)]);
                            
                            %Ensure that the alpha-peak location is included even under very low counts
                            indSel3 = unique([indSel3 find((E1 > mH*(pfit1d(4)*1000).^2/2/eV).*(E1 < 1.2*mA*(pfit1d(4)*1000).^2/2/2/eV) == 1)']);
                            
                            %Ensure that the limits stay consistent with the indices of the cropped box
                            indLim1 = [min(indSel1) max(indSel1)];
                            indLim2 = [min(indSel2) max(indSel2)];
                            indLim3 = [min(indSel3) max(indSel3)];

                            vHx = vHx(indSel1,indSel2,indSel3);      vAx = vAx(indSel1,indSel2,indSel3);
                            vHy = vHy(indSel1,indSel2,indSel3);      vAy = vAy(indSel1,indSel2,indSel3);
                            vHz = vHz(indSel1,indSel2,indSel3);      vAz = vAz(indSel1,indSel2,indSel3);
                            vHabs = vHabs(indSel1,indSel2,indSel3);  vAabs = vAabs(indSel1,indSel2,indSel3);
                
                            jq = jq(indSel1,indSel2,indSel3);
                            c = c(indSel1,indSel2,indSel3);
                            C = C(indSel1,indSel2,indSel3);
                            sigjq = sigjq(indSel1,indSel2,indSel3);
                            ph4 = ph4(indSel1,indSel2,indSel3);
                            th4 = th4(indSel1,indSel2,indSel3);
                            E4  =  E4(indSel1,indSel2,indSel3);
                            dE4 = dE4(indSel1,indSel2,indSel3);
                            Omega4 = Omega4(indSel1,indSel2,indSel3);
                            E1 = E1(indSel3);
                            c1 = c1(indSel3);

                        end
                    else
                        isCropped = false;
                    end

                    % Define phase-space sample locations for response function
                    maxdefv = 6.4; %Deflector/inner hemisphere voltage ratio for max. elevation
                    k_ideal = 7.8; %Ideal ESA k-factor

                    Vdef = maxdefv.*(squeeze(ph4(:,1,1))./deg2rad(45)); %Deflector voltages
                    
                    if indrSel > 0
    
                        fVdefIHall = SWIAr.fVdefIHall; 
                        fVdefIHall(SWIAr.def1or2all == 1) = -SWIAr.fVdefIHall(SWIAr.def1or2all == 1);
    
                        indResp = knnsearch(fVdefIHall',Vdef); %Indicies of the response matrix to use
                        
                        E4s = single(NaN([size(E4) size(SWIAr.RsAll,1)]));
                        ph4s = single(NaN([size(E4) size(SWIAr.RsAll,1)]));
                        th4s = single(NaN([size(E4) size(SWIAr.RsAll,1)]));
                        R4s  = single(NaN([size(E4) size(SWIAr.RsAll,1)]));
                        for indph = 1:length(Vdef)
                            indr = indResp(indph);
                            fE = SWIAr.ksAll(:,indr)./SWIAr.kAvg(indr); %Relative energies (Es/E)
                            dth = SWIAr.phvsAll(:,indr); %Relative azimuth angles
                            dph = SWIAr.thsAll(:,indr)-SWIAr.thAvg(indr); %Relative elevation angles phs-ph
                            
                            
                            %Sample locations in velocity space
                            E4s(indph,:,:,:) = repmat(E4(indph,:,:),[1 1 1 length(fE)]).*repmat(shiftdim(fE,-3),[size(E4(indph,:,:)) 1]);
                            th4s(indph,:,:,:) = repmat(th4(indph,:,:),[1 1 1 length(dth)]) + repmat(shiftdim(deg2rad(dth),-3),[size(E4(indph,:,:)) 1]);
                            ph4s(indph,:,:,:) = repmat(ph4(indph,:,:),[1 1 1 length(dph)]) + repmat(shiftdim(deg2rad(dph),-3),[size(E4(indph,:,:)) 1]);
                            R4s(indph,:,:,:) = repmat(shiftdim(SWIAr.RsAll(:,indr),-3),[size(E4(indph,:,:)) 1]);
                            
                        end
                    else
                        E4s = E4;
                        th4s = th4;
                        ph4s = ph4;
                        R4s = ones(size(E4));
                    end

                    vHsAbs = double(sqrt(2*E4s*eV/mH));
                    vAsAbs = double(sqrt(2*2*E4s*eV/mA));
                    [vHsx,vHsy,vHsz] = sph2cart(double(th4s),double(ph4s),vHsAbs);
                    [vAsx,vAsy,vAsz] = sph2cart(double(th4s),double(ph4s),vAsAbs);
                    R4s = double(R4s);
                                        
                    %Try fitting parameters for +B and -B beam 
                    for doingBplus = 0:1

                        %Skip either direction if the other is pre-selected
                        if (doingBplus == 0) && (SWIAf.useBeamBnegBpos(indf) == +1)
                            costfunMinBminus = inf;
                            continue;
                        elseif (doingBplus == 1) && (SWIAf.useBeamBnegBpos(indf) == -1)
                            costfunMinBplus = inf;
                            continue;
                        end

                        if doingBplus == 0
                            vnxSWIA = -double(SWIAf.BxnSWIA(indf));
                            vnySWIA = -double(SWIAf.BynSWIA(indf));
                            vnzSWIA = -double(SWIAf.BznSWIA(indf));
                        elseif doingBplus == 1
                            vnxSWIA = double(SWIAf.BxnSWIA(indf));
                            vnySWIA = double(SWIAf.BynSWIA(indf));
                            vnzSWIA = double(SWIAf.BznSWIA(indf));
                        else
                            error('Something went very wrong here');
                        end
                      
                      %Limit ub(13) beam diff. speed to the SWIA fine FOV
                      %Calculate the maximum beam speed within the SWIA fine FOV
                       Valf = Babs*1E-9/sqrt(mu0*(mH*(p0(1)+p0(3)) + mA*p0(2)));
                       vb = (0:1000);
                       vxb = p0(4) + vnxSWIA*vb;
                       vyb = p0(5) + vnySWIA*vb;
                       vzb = p0(6) + vnzSWIA*vb;
    
                       [azb,elb,Vb] = cart2sph(vxb,vyb,vzb);
                       azb(azb < 0) = azb(azb < 0)+2*pi;
    
                       indSelb = (azb > min(th4(:))).*(azb < max(th4(:))).*(elb > min(ph4(:))).*(elb < max(ph4(:)).*(Vb > min(vHabs(:)/1000)).*(Vb < max(vHabs(:)/1000))) == 1;
                       
                       if sum(indSelb) > 0
                           ub(13) = min([ub(13) max(vb(indSelb))/(Valf/1000)]);
                       end
                       if ub(13) < lb(13)
                           ub(13) = 1.1;
                       end
                    
                        %Estimate drift center of ring distribution
                        Vx = pfit1d(4)*(SWIAf.VxSWIA(indf)/SWIAf.Vabs(indf));
                        Vy = pfit1d(4)*(SWIAf.VySWIA(indf)/SWIAf.Vabs(indf));
                        Vz = pfit1d(4)*(SWIAf.VzSWIA(indf)/SWIAf.Vabs(indf));
                        
                        Vrx = double(Vx - Bxn*(Vx*Bxn + Vy*Byn + Vz*Bzn)/(Bxn^2 + Byn^2 + Bzn^2));
                        Vry = double(Vy - Byn*(Vx*Bxn + Vy*Byn + Vz*Bzn)/(Bxn^2 + Byn^2 + Bzn^2));
                        Vrz = double(Vz - Bzn*(Vx*Bxn + Vy*Byn + Vz*Bzn)/(Bxn^2 + Byn^2 + Bzn^2));
                        
                        VrAbs = sqrt(Vrx^2 + Vry^2 + Vrz^2);
    
                        % Parallel and perpendicular velocity components
                        funvHpara = @(p) ((vHsx-p(4)*1000).*Bxn + (vHsy-p(5)*1000).*Byn + (vHsz-p(6)*1000).*Bzn);
                        funvHperp = @(p) sqrt(((vHsx-p(4)*1000)-funvHpara(p).*Bxn).^2 + ((vHsy-p(5)*1000)-funvHpara(p).*Byn).^2 + ((vHsz-p(6)*1000)-funvHpara(p).*Bzn).^2);
                        
                        funvApara = @(p) ((vAsx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn).*Bxn + (vAsy-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn).*Byn + (vAsz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn).*Bzn);
                        funvAperp = @(p) sqrt(((vAsx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn)-funvApara(p).*Bxn).^2 + ((vAsy-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn)-funvApara(p).*Byn).^2 + ((vAsz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn)-funvApara(p).*Bzn).^2);
                        
                        funvHparaRing = @(p) ((vHsx-Vrx*1000).*Bxn + (vHsy-Vry*1000).*Byn + (vHsz-Vrz*1000).*Bzn);
                        funvHperpRing = @(p) sqrt(((vHsx-Vrx*1000)-funvHparaRing(p).*Bxn).^2 + ((vHsy-Vry*1000)-funvHparaRing(p).*Byn).^2 + ((vHsz-Vrz*1000)-funvHparaRing(p).*Bzn).^2);
                        
                        % Distribution function
                        if indrSel >= (length(SWIArc)-1)
                            Tring = 3E+4; %Ring temperature [K]
                        else
                            Tring = 1E+5; %Artificially large if not using high-resolution response function
                        end
    
                        funA = @(p) exp(-(VrAbs*1000)^2/sqrt(2*k*Tring/mH)^2) + sqrt(pi)*((VrAbs*1000)/sqrt(2*k*Tring/mH))*erfc(-(VrAbs*1000)/sqrt(2*k*Tring/mH)); %Normalization parameter
                        jfuns  = @(p) 1E-4*vHsAbs.^2./mH.*eV.*p(1).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(10))/(3*p(10)^(2/3)))^(3/2)*sqrt(2*k*p(7)/mH)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(10))/(3*sqrt(2*k*p(7)/mH)^2)*(funvHpara(p).^2 + 1/p(10)*funvHperp(p).^2)).^(-p(12)-1)    + 1E-4*vHsAbs.^2./mH.*eV.*p(3).*(pi*sqrt(2*k*p(7)/mH)^2).^(-3/2).*(p(12)-3/2).^(-3/2).*gamma(p(12)+1)./gamma(p(12)-1/2).*(1 + ((vHsx-(p(4)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnxSWIA)).^2 + (vHsy-(p(5)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnySWIA)).^2 + (vHsz-(p(6)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnzSWIA)).^2)./((p(12)-3/2).*sqrt(2*k*p(7)/mH).^2)).^(-p(12)-1)     + 2*1E-4*vAsAbs.^2./mA.*eV.*p(2).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(11))/(3*p(11)^(2/3)))^(3/2)*sqrt(2*k*p(9)*p(7)/mA)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(11))/(3*sqrt(2*k*p(9)*p(7)/mA)^2)*(funvApara(p).^2 + 1/p(11)*funvAperp(p).^2)).^(-p(12)-1)    + 1E-4*vHsAbs.^2./mH.*eV*p(14)*(pi^(3/2)*sqrt(2*k*Tring/mH)^3*funA(p))^(-1).*exp(-(funvHparaRing(p).^2)/(sqrt(2*k*Tring/mH))^2).*exp(-(funvHperpRing(p)-(VrAbs*1000)).^2/(sqrt(2*k*Tring/mH))^2);
                        jfun = @(p) sum(jfuns(p).*R4s,4);
    
                        % Cost function based on log-Poisson PMF and measured counts 
                        costfun  = @(p) -sum(logPoisspdf(c,jfun(p).*E4./C+p(15)),'all'); %Probability cost function 
                        
                        % Optimization
                        opts = optimoptions('fmincon','Display','off','Algorithm','interior-point','MaxFunctionEvaluations',8000,'OptimalityTolerance',1E-11);
                        [pfit,costfunMin,exitFlag,~,~,~,H] = fmincon(costfun,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
                        
                        if doingBplus == 0
                            vnSWIABminus = [vnxSWIA,vnySWIA,vnzSWIA];
                            pfitBminus = pfit;
                            costfunMinBminus = costfunMin;
                            exitFlagBminus = exitFlag;
                            HBminus = H;

                            jfunsBminus = jfuns;
                            jfunBminus = jfun;

                        elseif doingBplus == 1
                            vnSWIABplus = [vnxSWIA,vnySWIA,vnzSWIA];
                            pfitBplus = pfit;
                            costfunMinBplus = costfunMin;
                            exitFlagBplus = exitFlag;
                            HBplus = H;

                            jfunsBplus = jfuns;
                            jfunBplus = jfun;
                        end
                        
                    end
                    
                    %Select the most probable fit
                    if costfunMinBminus < costfunMinBplus
                        vnxSWIA = vnSWIABminus(1);
                        vnySWIA = vnSWIABminus(2);
                        vnzSWIA = vnSWIABminus(3);
                        pfit = pfitBminus;
                        costfunMin = costfunMinBminus;
                        exitFlag = exitFlagBminus;
                        H = HBminus;

                        jfuns = jfunsBminus;
                        jfun = jfunBminus;
                    else
                        vnxSWIA = vnSWIABplus(1);
                        vnySWIA = vnSWIABplus(2);
                        vnzSWIA = vnSWIABplus(3);
                        pfit = pfitBplus;
                        costfunMin = costfunMinBplus;
                        exitFlag = exitFlagBplus;
                        H = HBplus;
                        
                        jfuns = jfunsBplus;
                        jfun = jfunBplus;
                    end

                    I = H; % Fisher information matrix (= Hessian when using -log MLE)
                    cov = inv(I); % Lower limit for the covariance matrix (Cramer-Rao bound)
                    sigpfit = sqrt(diag(cov)); % Uncertainties of fitted parameters

                    % Calculate bulk parameters
                    NH = (pfit(1)+pfit(3))*1E-6; %Total solar wind proton densitiy (Core+Beam)
                    NA =  pfit(2)*1E-6; %Solar wind He++ density
                    NHring = pfit(13);
                    
                    NHcore = pfit(1)*1E-6;
                    VHcorexSWIA = pfit(4);
                    VHcoreySWIA = pfit(5);
                    VHcorezSWIA = pfit(6);

                    %Beam velocities
                    Valf = Babs*1E-9/sqrt(mu0*(mH*(pfit(1) + pfit(3)) + mA*pfit(2))); %Alfv?n wave speed
                    
                    dVbeamxSWIA = 1E-3*pfit(13)*Valf*vnxSWIA; %Differential beam velocity (relative core)
                    dVbeamySWIA = 1E-3*pfit(13)*Valf*vnySWIA;
                    dVbeamzSWIA = 1E-3*pfit(13)*Valf*vnzSWIA;

                    NHbeam = pfit(3)*1E-6;
                    VHbeamxSWIA = VHcorexSWIA + dVbeamxSWIA;
                    VHbeamySWIA = VHcoreySWIA + dVbeamySWIA;
                    VHbeamzSWIA = VHcorezSWIA + dVbeamzSWIA;
                    
                    %Net H+ velocity as density-weighted velocities of core and beam
                    VHxSWIA = (NHcore.*VHcorexSWIA + NHbeam*VHbeamxSWIA)/(NHcore + NHbeam); 
                    VHySWIA = (NHcore.*VHcoreySWIA + NHbeam*VHbeamySWIA)/(NHcore + NHbeam); 
                    VHzSWIA = (NHcore.*VHcorezSWIA + NHbeam*VHbeamzSWIA)/(NHcore + NHbeam);
                    
                    VAxSWIA = pfit(4) + 1E-3*(pfit(8)*Valf*Bxn);
                    VAySWIA = pfit(5) + 1E-3*(pfit(8)*Valf*Byn);
                    VAzSWIA = pfit(6) + 1E-3*(pfit(8)*Valf*Bzn);
                    
                    TH = pfit(7); %Proton scalar temperature
                    TA = pfit(9)*pfit(7); %He++ scalar temperature
                    alphaH = pfit(10); %Alpha parameter (temperature anisotropy factor)
                    alphaA = pfit(11); %Alpha parameter (temperature anisotropy factor)
                    kappa = pfit(12); %Kappa parameter (thermalization factor)
                    
                    %If very low kappa, redo with highest possible resolution response matrix
                    if ((kappa < 1.8) || (TH < 8E+3) || (TA < 1E+4)) && (indrSel < length(SWIArc)) && isnan(SWIAf.useRespMat(indf))
                        [~,elVHSWIA,~] = cart2sph(VHxSWIA,VHySWIA,VHzSWIA);
                        if abs(rad2deg(elVHSWIA)) > 30
                            SWIAf.useRespMat(indf) = length(SWIArc); %Widely sampled at high elevations, but oversampled in center
                        else
                            SWIAf.useRespMat(indf) = length(SWIArc)-1;
                        end
                        %Re-do the distribution if the forced response matrix has a higher resolution
                        if indrSel < SWIAf.useRespMat(indf)
                            jj = jj-1;
                            indsFineProc(indsFineProc == indf) = []; %De-log that this distribution has been processed
                            continue;
                        else
                            SWIAf.useRespMat(indf) = NaN;
                        end
                    end
                    
                    %Azimuth, elevation, and speed of the bulk velocity vectors
                    [azHswia,elHswia,VHabs] = cart2sph(VHxSWIA,VHySWIA,VHzSWIA);
                    [azAswia,elAswia,VAabs] = cart2sph(VAxSWIA,VAySWIA,VAzSWIA);
                    
                    %Convert to MSO
                    tempH = SWIAf.Mswia2msoAll(:,:,indf)*[VHxSWIA;VHySWIA;VHzSWIA];
                    VHxMSO = tempH(1); VHyMSO = tempH(2); VHzMSO = tempH(3);
                    
                    tempHcore = SWIAf.Mswia2msoAll(:,:,indf)*[VHcorexSWIA;VHcoreySWIA;VHcorezSWIA];
                    VHcorexMSO = tempHcore(1); VHcoreyMSO = tempHcore(2); VHcorezMSO = tempHcore(3);

                    tempHdbeam = SWIAf.Mswia2msoAll(:,:,indf)*[dVbeamxSWIA;dVbeamySWIA;dVbeamzSWIA];
                    dVHbeamxMSO = tempHdbeam(1); dVHbeamyMSO = tempHdbeam(2); dVHbeamzMSO = tempHdbeam(3);
                    
                    tempA = SWIAf.Mswia2msoAll(:,:,indf)*[VAxSWIA;VAySWIA;VAzSWIA];
                    VAxMSO = tempA(1); VAyMSO = tempA(2); VAzMSO = tempA(3);
                    
                    %Uncertainties

                    %Identify edge solutions
                    if sum(imag(sigpfit) > 0)
                        warning('Bad Hessian from 3D fit, solution likely at boundary');
                        clearvars sigpfit;
                        sigpfit = inf(14,1);
                    end

                    sigNH = sqrt(sigpfit(1)^2+sigpfit(3)^2+2*cov(1,3))*1E-6;
                    sigNHcore = sqrt(sigpfit(1)^2);
                    sigNHbeam = sqrt(sigpfit(3)^2);
                    sigNA = sigpfit(2)*1E-6;

                    sigVHcorexSWIA = sqrt(sigpfit(4)^2);
                    sigVHcoreySWIA = sqrt(sigpfit(5)^2);
                    sigVHcorezSWIA = sqrt(sigpfit(6)^2);

                    sigdVHbeamxSWIA = sqrt((1E-3*Valf*vnxSWIA)^2*sigpfit(13)^2);
                    sigdVHbeamySWIA = sqrt((1E-3*Valf*vnySWIA)^2*sigpfit(13)^2);
                    sigdVHbeamzSWIA = sqrt((1E-3*Valf*vnzSWIA)^2*sigpfit(13)^2);

                    sigVHxSWIA = sqrt((NHcore/NH)^2*sigVHcorexSWIA^2 + (NHbeam/NH)^2*sigdVHbeamxSWIA^2);
                    sigVHySWIA = sqrt((NHcore/NH)^2*sigVHcoreySWIA^2 + (NHbeam/NH)^2*sigdVHbeamySWIA^2);
                    sigVHzSWIA = sqrt((NHcore/NH)^2*sigVHcorezSWIA^2 + (NHbeam/NH)^2*sigdVHbeamzSWIA^2);
                    
                    sigVAxSWIA = sqrt(sigVHcorexSWIA^2 + (Valf*Bxn)^2*sigpfit(8)^2);% + (pfit(8)*Bxn)^2*sigValf^2);
                    sigVAySWIA = sqrt(sigVHcoreySWIA^2 + (Valf*Byn)^2*sigpfit(8)^2);
                    sigVAzSWIA = sqrt(sigVHcorezSWIA^2 + (Valf*Bzn)^2*sigpfit(8)^2);

                    sigVHcorexMSO = sqrt(SWIAf.Mswia2msoAll(1,1,indf).^2*sigVHcorexSWIA.^2 + SWIAf.Mswia2msoAll(1,2,indf).^2*sigVHcoreySWIA.^2 + SWIAf.Mswia2msoAll(1,3,indf).^2*sigVHcorezSWIA.^2);
                    sigVHcoreyMSO = sqrt(SWIAf.Mswia2msoAll(2,1,indf).^2*sigVHcorexSWIA.^2 + SWIAf.Mswia2msoAll(2,2,indf).^2*sigVHcoreySWIA.^2 + SWIAf.Mswia2msoAll(2,3,indf).^2*sigVHcorezSWIA.^2);
                    sigVHcorezMSO = sqrt(SWIAf.Mswia2msoAll(3,1,indf).^2*sigVHcorexSWIA.^2 + SWIAf.Mswia2msoAll(3,2,indf).^2*sigVHcoreySWIA.^2 + SWIAf.Mswia2msoAll(3,3,indf).^2*sigVHcorezSWIA.^2);
                    
                    sigdVHbeamxMSO = sqrt(SWIAf.Mswia2msoAll(1,1,indf).^2*sigdVHbeamxSWIA.^2 + SWIAf.Mswia2msoAll(1,2,indf).^2*sigdVHbeamySWIA.^2 + SWIAf.Mswia2msoAll(1,3,indf).^2*sigdVHbeamzSWIA.^2);
                    sigdVHbeamyMSO = sqrt(SWIAf.Mswia2msoAll(2,1,indf).^2*sigdVHbeamxSWIA.^2 + SWIAf.Mswia2msoAll(2,2,indf).^2*sigdVHbeamySWIA.^2 + SWIAf.Mswia2msoAll(2,3,indf).^2*sigdVHbeamzSWIA.^2);
                    sigdVHbeamzMSO = sqrt(SWIAf.Mswia2msoAll(3,1,indf).^2*sigdVHbeamxSWIA.^2 + SWIAf.Mswia2msoAll(3,2,indf).^2*sigdVHbeamySWIA.^2 + SWIAf.Mswia2msoAll(3,3,indf).^2*sigdVHbeamzSWIA.^2);
                    
                    sigVHxMSO = sqrt(SWIAf.Mswia2msoAll(1,1,indf).^2*sigVHxSWIA.^2 + SWIAf.Mswia2msoAll(1,2,indf).^2*sigVHySWIA.^2 + SWIAf.Mswia2msoAll(1,3,indf).^2*sigVHzSWIA.^2);
                    sigVHyMSO = sqrt(SWIAf.Mswia2msoAll(2,1,indf).^2*sigVHxSWIA.^2 + SWIAf.Mswia2msoAll(2,2,indf).^2*sigVHySWIA.^2 + SWIAf.Mswia2msoAll(2,3,indf).^2*sigVHzSWIA.^2);
                    sigVHzMSO = sqrt(SWIAf.Mswia2msoAll(3,1,indf).^2*sigVHxSWIA.^2 + SWIAf.Mswia2msoAll(3,2,indf).^2*sigVHySWIA.^2 + SWIAf.Mswia2msoAll(3,3,indf).^2*sigVHzSWIA.^2);
                    
                    sigVAxMSO = sqrt(SWIAf.Mswia2msoAll(1,1,indf).^2*sigVAxSWIA.^2 + SWIAf.Mswia2msoAll(1,2,indf).^2*sigVAySWIA.^2 + SWIAf.Mswia2msoAll(1,3,indf).^2*sigVAzSWIA.^2);
                    sigVAyMSO = sqrt(SWIAf.Mswia2msoAll(2,1,indf).^2*sigVAxSWIA.^2 + SWIAf.Mswia2msoAll(2,2,indf).^2*sigVAySWIA.^2 + SWIAf.Mswia2msoAll(2,3,indf).^2*sigVAzSWIA.^2);
                    sigVAzMSO = sqrt(SWIAf.Mswia2msoAll(3,1,indf).^2*sigVAxSWIA.^2 + SWIAf.Mswia2msoAll(3,2,indf).^2*sigVAySWIA.^2 + SWIAf.Mswia2msoAll(3,3,indf).^2*sigVAzSWIA.^2);
              
                    sigTH = sigpfit(7);
                    sigTA = sqrt(TA^2*((sigpfit(9)/pfit(9))^2 + (sigpfit(7)/pfit(7))^2 + 2*cov(9,7)/(pfit(7)*pfit(9))));
                    sigAlphaH = sigpfit(10);
                    sigAlphaA = sigpfit(11);
                    sigKappa = sigpfit(12);
                    
                    %Pad distribution to 48 E-indices if SVY or otherwise cropped
                    if SWIAf.grouping(indf) == 0 && ~isCropped
                        Jq = squeeze(sum(sum(jfun(pfit).*Omega4,2),1));
                    elseif SWIAf.grouping(indf) ~= 0 && ~isCropped
                        Jq = [NaN(8,1);squeeze(sum(sum(jfun(pfit).*Omega4,2),1));NaN(8,1)];
                    elseif SWIAf.grouping(indf) == 0 && isCropped
                        Jq = [NaN(indLim3(1)-1,1);squeeze(sum(sum(jfun(pfit).*Omega4,2),1));NaN(48-indLim3(2),1)];
                    end

                    if length(Jq) ~= 48
                        disp(['   %% Warning: ii=',num2str(ii),': length(Jq) ~= 48, skipping']);
                        continue;
                    end
                    
                    %Estimate goodness of fit (Normal approximation)
                    cm4 = jfun(pfit).*E4./C+pfit(15);
                    frac95 = sum(sum(sum((c >= poissinv(1-0.95,cm4)).*(c <= poissinv(0.95,cm4)))))/numel(c);
                    
                    ind = cm4 >= 10; 
                    chi2rApprox = sum((cm4(ind)-c(ind)).^2./sqrt(cm4(ind)+0.5).^2,'all')/numel(c(ind));

                    %Analyze whether the distribution is likely undisturbed solar wind
                    isSW = true;
                    VthH = sqrt(2*k*TH/mH); 
                    
                    %Thermal speed to high (likely sheath)
                    if VthH/(1000*VHabs) > 0.2
                        isSW = false;
                    end
                    
                    %Exclude if too deflected
                    VHxnSWIA = VHxSWIA/VHabs;
                    VHynSWIA = VHySWIA/VHabs;
                    VHznSWIA = VHzSWIA/VHabs;
                    dthvsw = acosd(SWIAf.vswnSWIA(1,indf)*VHxnSWIA + SWIAf.vswnSWIA(2,indf)*VHynSWIA + SWIAf.vswnSWIA(3,indf)*VHznSWIA);
                    if dthvsw > 10
                        isSW = false;
                    end

                    %Go to Coarse fits if bulk outside or at edge of Fine FOV (sqrt(4.5^2 + 3.75^2) = 5.86 deg)
                    [vxn,vyn,vzn] = sph2cart(th4(2:(end-1),2:(end-1),round(size(c,3)/2)),ph4(2:(end-1),2:(end-1),round(size(c,3)/2)),1);
                    dthFOV = acosd(vxn.*VHxnSWIA + vyn.*VHynSWIA + vzn.*VHznSWIA);
                    if min(dthFOV(:)) > 6
                        indsFineIncl(jj) = false;
                        continue;
                    end

                    %Override if in either inclusion or exclusion interval
                    if SWIAf.indIncl(indf)
                        isSW = true;
                    end
                    if SWIAf.indExcl(indf)
                        isSW = false;
                    end

                    %% Generate deconvoluted (natural) spectra for reference
                        dazel = 0.5; %Angular resolution 
                        En = 10.^(linspace(log10(SWIAf.E(end,indf)),log10(SWIAf.E(1,indf)),numelsE));
                        az1 = 130:dazel:230;
                        el1 = -55:dazel:55;
                        [azg,elg,Eqg] = ndgrid(az1,el1,En);
                        [~,~,dEqp] = ndgrid(az1,el1,gradient(En));
                        [vHnx,vHny,vHnz] = sph2cart(deg2rad(azg),deg2rad(elg),sqrt(2*Eqg*eV/mH));
                        vHnAbs = sqrt(vHnx.^2 + vHny.^2 + vHnz.^2);
                        
                        [vAnx,vAny,vAnz] = sph2cart(deg2rad(azg),deg2rad(elg),sqrt(2*2*Eqg*eV/mA));
                        vAnAbs = sqrt(vAnx.^2 + vAny.^2 + vAnz.^2);

                        clearvars funnvHpara funnvHperp funnvApara funnvAperp;
                        funnvHpara = @(p) ((vHnx-p(4)*1000).*Bxn + (vHny-p(5)*1000).*Byn + (vHnz-p(6)*1000).*Bzn);
                        funnvHperp = @(p) sqrt(((vHnx-p(4)*1000)-funnvHpara(p).*Bxn).^2 + ((vHny-p(5)*1000)-funnvHpara(p).*Byn).^2 + ((vHnz-p(6)*1000)-funnvHpara(p).*Bzn).^2);
                        
                        funnvApara = @(p) ((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn).*Bxn + (vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn).*Byn + (vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn).*Bzn);
                        funnvAperp = @(p) sqrt(((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn)-funnvApara(p).*Bxn).^2 + ((vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn)-funnvApara(p).*Byn).^2 + ((vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn)-funnvApara(p).*Bzn).^2);
                        
                        funnvHparaRing = @(p) ((vHnx-Vrx*1000).*Bxn + (vHny-Vry*1000).*Byn + (vHnz-Vrz*1000).*Bzn);
                        funnvHperpRing = @(p) sqrt(((vHnx-Vrx*1000)-funnvHparaRing(p).*Bxn).^2 + ((vHny-Vry*1000)-funnvHparaRing(p).*Byn).^2 + ((vHnz-Vrz*1000)-funnvHparaRing(p).*Bzn).^2);
                    
                        funA = @(p) exp(-(VrAbs*1000)^2/sqrt(2*k*Tring/mH)^2) + sqrt(pi)*((VrAbs*1000)/sqrt(2*k*Tring/mH))*erfc(-(VrAbs*1000)/sqrt(2*k*Tring/mH)); %Normalization parameter
                        
                        jfunnCore = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(1).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(10))/(3*p(10)^(2/3)))^(3/2)*sqrt(2*k*p(7)/mH)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(10))/(3*sqrt(2*k*p(7)/mH)^2)*(funnvHpara(p).^2 + 1/p(10)*funnvHperp(p).^2)).^(-p(12)-1);
                        jfunnBeam = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(3).*(pi*sqrt(2*k*p(7)/mH)^2).^(-3/2).*(p(12)-3/2).^(-3/2).*gamma(p(12)+1)./gamma(p(12)-1/2).*(1 + ((vHnx-(p(4)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnxSWIA)).^2 + (vHny-(p(5)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnySWIA)).^2 + (vHnz-(p(6)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnzSWIA)).^2)./((p(12)-3/2).*sqrt(2*k*p(7)/mH).^2)).^(-p(12)-1);
                        jfunnA    = @(p) 2*1E-4*vAnAbs.^2./mA.*eV.*p(2).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(11))/(3*p(11)^(2/3)))^(3/2)*sqrt(2*k*p(9)*p(7)/mA)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(11))/(3*sqrt(2*k*p(9)*p(7)/mA)^2)*(funnvApara(p).^2 + 1/p(11)*funnvAperp(p).^2)).^(-p(12)-1);
                        jfunnRing = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(14)*(pi^(3/2)*sqrt(2*k*Tring/mH)^3*funA(p))^(-1).*exp(-(funnvHparaRing(p).^2)/(sqrt(2*k*Tring/mH))^2).*exp(-(funnvHperpRing(p)-(VrAbs*1000)).^2/(sqrt(2*k*Tring/mH))^2);
                                                
                        Omegang = deg2rad(dazel).^2.*cosd(elg);
                        jEnCoreOmni = squeeze(sum(sum(jfunnCore(pfit).*Eqg.*Omegang,2),1))'; 
                        jEnBeamOmni = squeeze(sum(sum(jfunnBeam(pfit).*Eqg.*Omegang,2),1))';
                        jEnAomni = squeeze(sum(sum(jfunnA(pfit).*Eqg.*Omegang,2),1))';
                        jEnRingOmni = squeeze(sum(sum(jfunnRing(pfit).*Eqg.*Omegang,2),1))';  
                        
                        %Estimate bulk (core+beam) deconvoluted proton temperature by integrating the full H+ distribution
                        jHEnTot = jfunnCore(pfit)+jfunnBeam(pfit); %Total solar wind deconvoluted H differential flux

                        Pxx = mH*squeeze(nansum(nansum(nansum(jHEnTot.*10^4.*vHnAbs.*cosd(azg).^2.*cosd(elg).^2.*Omegang.*dEqp,3),2),1)) - mH*NH.*1E+6.*(VHxSWIA*1000).^2;
                        Pyy = mH*squeeze(nansum(nansum(nansum(jHEnTot.*10^4.*vHnAbs.*sind(azg).^2.*cosd(elg).^2.*Omegang.*dEqp,3),2),1)) - mH*NH.*1E+6.*(VHySWIA*1000).^2;
                        Pzz = mH*squeeze(nansum(nansum(nansum(jHEnTot.*10^4.*vHnAbs.*sind(elg).^2.*Omegang.*dEqp,3),2),1)) - mH*NH.*1E+6.*(VHzSWIA*1000).^2;
                        
                        %Derivation of parallel/perpendicular temperatures
                        %following Verscharen 2019 (eqs. 36-39) doi:10.1007/s41116-019-0021-0

                        Ppara = Pxx*Bxn.^2 + Pyy*Byn.^2 + Pzz*Bzn.^2;
                        Pperp = Pxx*(1-Bxn^2)/2 + Pyy*(1-Byn^2)/2 + Pzz*(1-Bzn^2)/2;
                        
                        THtotPara = Ppara/(NH*1E+6*k);
                        THtotPerp = Pperp/(NH*1E+6*k);
                        THtot = abs(2*Pperp + Ppara)./(3*NH*1E+6*k);


%%
                    % Concatenate results of fit into common structure
                    SWIA.time      = cat(1,SWIA.time,SWIAf.time(indf));
                    SWIA.NH        = cat(1,SWIA.NH,NH);
                    SWIA.NHcore    = cat(1,SWIA.NHcore,NHcore);
                    SWIA.NHbeam    = cat(1,SWIA.NHbeam,NHbeam);
                    SWIA.NA        = cat(1,SWIA.NA,NA);
                    SWIA.VHx       = cat(1,SWIA.VHx,VHxMSO);
                    SWIA.VHy       = cat(1,SWIA.VHy,VHyMSO);
                    SWIA.VHz       = cat(1,SWIA.VHz,VHzMSO);
                    SWIA.VHabs     = cat(1,SWIA.VHabs,VHabs);
                    SWIA.VAx       = cat(1,SWIA.VAx,VAxMSO);
                    SWIA.VAy       = cat(1,SWIA.VAy,VAyMSO);
                    SWIA.VAz       = cat(1,SWIA.VAz,VAzMSO);
                    SWIA.VAabs     = cat(1,SWIA.VAabs,VAabs);
                    SWIA.THcore    = cat(1,SWIA.THcore,TH);
                    SWIA.TA        = cat(1,SWIA.TA,TA);
                    SWIA.TH        = cat(1,SWIA.TH,THtot);
                    SWIA.THperp    = cat(1,SWIA.THperp,THtotPerp);
                    SWIA.THpara    = cat(1,SWIA.THpara,THtotPara);
                    SWIA.alphaHcore = cat(1,SWIA.alphaHcore,alphaH);
                    SWIA.alphaA     = cat(1,SWIA.alphaA,alphaA);
                    SWIA.kappa     = cat(1,SWIA.kappa,kappa);

                    SWIA.VHcorex       = cat(1,SWIA.VHcorex,VHcorexMSO);
                    SWIA.VHcorey       = cat(1,SWIA.VHcorey,VHcoreyMSO);
                    SWIA.VHcorez       = cat(1,SWIA.VHcorez,VHcorezMSO);

                    SWIA.dVHbeamx       = cat(1,SWIA.dVHbeamx,dVHbeamxMSO);
                    SWIA.dVHbeamy       = cat(1,SWIA.dVHbeamy,dVHbeamyMSO);
                    SWIA.dVHbeamz       = cat(1,SWIA.dVHbeamz,dVHbeamzMSO);

                    SWIA.vnxSWIA   = cat(1,SWIA.vnxSWIA,vnxSWIA);
                    SWIA.vnySWIA   = cat(1,SWIA.vnySWIA,vnySWIA);
                    SWIA.vnzSWIA   = cat(1,SWIA.vnzSWIA,vnzSWIA);

                    SWIA.Bx        = cat(1,SWIA.Bx,SWIAf.Bx(indf));
                    SWIA.By        = cat(1,SWIA.By,SWIAf.By(indf));
                    SWIA.Bz        = cat(1,SWIA.Bz,SWIAf.Bz(indf));
                    
                    SWIA.BxSWIA    = cat(1,SWIA.BxSWIA,SWIAf.BxSWIA(indf));
                    SWIA.BySWIA    = cat(1,SWIA.BySWIA,SWIAf.BySWIA(indf));
                    SWIA.BzSWIA    = cat(1,SWIA.BzSWIA,SWIAf.BzSWIA(indf));

                    SWIA.Valf      = cat(1,SWIA.Valf,Valf);

                    SWIA.En           = cat(1,SWIA.En,En);
                    SWIA.jEnCoreOmni = cat(1,SWIA.jEnCoreOmni,jEnCoreOmni);
                    SWIA.jEnBeamOmni = cat(1,SWIA.jEnBeamOmni,jEnBeamOmni);
                    SWIA.jEnAomni = cat(1,SWIA.jEnAomni,jEnAomni);
                    SWIA.jEnRingOmni = cat(1,SWIA.jEnRingOmni,jEnRingOmni);

                    SWIA.sigNH     = cat(1,SWIA.sigNH,sigNH);
                    SWIA.sigNA     = cat(1,SWIA.sigNA,sigNA);

                    SWIA.sigNHcore     = cat(1,SWIA.sigNHcore,sigNHcore);
                    SWIA.sigNHbeam     = cat(1,SWIA.sigNHbeam,sigNHbeam);
                    
                    SWIA.sigVHcorex = cat(1,SWIA.sigVHcorex,sigVHcorexMSO);
                    SWIA.sigVHcorey = cat(1,SWIA.sigVHcorey,sigVHcoreyMSO);
                    SWIA.sigVHcorez = cat(1,SWIA.sigVHcorez,sigVHcorezMSO);

                    SWIA.sigdVHbeamx = cat(1,SWIA.sigdVHbeamx,sigdVHbeamxMSO);
                    SWIA.sigdVHbeamy = cat(1,SWIA.sigdVHbeamy,sigdVHbeamyMSO);
                    SWIA.sigdVHbeamz = cat(1,SWIA.sigdVHbeamz,sigdVHbeamzMSO);

                    SWIA.sigVHx = cat(1,SWIA.sigVHx,sigVHxMSO);
                    SWIA.sigVHy = cat(1,SWIA.sigVHy,sigVHyMSO);
                    SWIA.sigVHz = cat(1,SWIA.sigVHz,sigVHzMSO);
                    SWIA.sigVAx = cat(1,SWIA.sigVAx,sigVAxMSO);
                    SWIA.sigVAy = cat(1,SWIA.sigVAy,sigVAyMSO);
                    SWIA.sigVAz = cat(1,SWIA.sigVAz,sigVAzMSO);

                    SWIA.sigTHcore     = cat(1,SWIA.sigTHcore,sigTH);
                    SWIA.sigTA         = cat(1,SWIA.sigTA,sigTA);
                    SWIA.sigAlphaHcore = cat(1,SWIA.sigAlphaHcore,sigAlphaH);
                    SWIA.sigAlphaA     = cat(1,SWIA.sigAlphaA,sigAlphaA);
                    SWIA.sigKappa      = cat(1,SWIA.sigKappa,sigKappa);

                    SWIA.Eq          = cat(1,SWIA.Eq,SWIAf.E(:,indf)');
                    SWIA.Jq          = cat(1,SWIA.Jq,Jq');
                    SWIA.pfit        = cat(1,SWIA.pfit,pfit);
                    SWIA.indrSel     = cat(1,SWIA.indrSel,indrSel);
                    SWIA.atten_state = cat(1,SWIA.atten_state,SWIAf.atten_state(indf));
                    SWIA.frac95      = cat(1,SWIA.frac95,frac95);
                    SWIA.chi2rApprox = cat(1,SWIA.chi2rApprox,chi2rApprox);

                    SWIA.Mswia2mso = cat(1,SWIA.Mswia2mso,shiftdim(SWIAf.Mswia2msoAll(:,:,indf),-1));
                    SWIA.Mmso2swia = cat(1,SWIA.Mmso2swia,shiftdim(SWIAf.Mmso2swiaAll(:,:,indf),-1));
                    SWIA.isSW      = cat(1,SWIA.isSW,isSW);
                    SWIA.isDraped  = cat(1,SWIA.isDraped,SWIAc.isDraped(ii));
                    SWIA.isIMF     = cat(1,SWIA.isIMF,SWIAc.isIMF(ii));
                    SWIA.isIono    = cat(1,SWIA.isIono,SWIAc.isIono(ii));
                    SWIA.isFine    = cat(1,SWIA.isFine,true);
                    SWIA.isFS      = cat(1,SWIA.isFS,SWIAc.isFS(ii));
                    SWIA.isSvy     = cat(1,SWIA.isSvy,SWIAf.grouping(indf) == 1);
                    SWIA.fitExitFlag  = cat(1,SWIA.fitExitFlag,exitFlag);
                    
                    SWIA.xSC     = cat(1,SWIA.xSC,SWIAf.xSC(indf));
                    SWIA.ySC     = cat(1,SWIA.ySC,SWIAf.ySC(indf));
                    SWIA.zSC     = cat(1,SWIA.zSC,SWIAf.zSC(indf));

                    SWIA.dtRunTime = cat(1,SWIA.dtRunTime,toc);

                    indsFineLog = cat(1,indsFineLog,indf);

                    if doplot == 1
                        %%
                        jfunsCore = @(p)   1E-4*vHsAbs.^2./mH.*eV.*p(1).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(10))/(3*p(10)^(2/3)))^(3/2)*sqrt(2*k*p(7)/mH)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(10))/(3*sqrt(2*k*p(7)/mH)^2)*(funvHpara(p).^2 + 1/p(10)*funvHperp(p).^2)).^(-p(12)-1);
                        jfunsBeam = @(p)   1E-4*vHsAbs.^2./mH.*eV.*p(3).*(pi*sqrt(2*k*p(7)/mH)^2).^(-3/2).*(p(12)-3/2).^(-3/2).*gamma(p(12)+1)./gamma(p(12)-1/2).*(1 + ((vHsx-(p(4)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnxSWIA)).^2 + (vHsy-(p(5)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnySWIA)).^2 + (vHsz-(p(6)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnzSWIA)).^2)./((p(12)-3/2).*sqrt(2*k*p(7)/mH).^2)).^(-p(12)-1);
                        jfunsA    = @(p) 2*1E-4*vAsAbs.^2./mA.*eV.*p(2).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(11))/(3*p(11)^(2/3)))^(3/2)*sqrt(2*k*p(9)*p(7)/mA)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(11))/(3*sqrt(2*k*p(9)*p(7)/mA)^2)*(funvApara(p).^2 + 1/p(11)*funvAperp(p).^2)).^(-p(12)-1);
                        jfunsRing = @(p)   1E-4*vHsAbs.^2./mH.*eV*p(14)*(pi^(3/2)*sqrt(2*k*Tring/mH)^3*funA(p))^(-1).*exp(-(funvHparaRing(p).^2)/(sqrt(2*k*Tring/mH))^2).*exp(-(funvHperpRing(p)-(VrAbs*1000)).^2/(sqrt(2*k*Tring/mH))^2);
                                          
                        jfunCore = @(p) sum(jfunsCore(p).*R4s,4);
                        jfunBeam = @(p) sum(jfunsBeam(p).*R4s,4);
                           jfunA = @(p) sum(jfunsA(p).*R4s,4);
                        jfunRing = @(p) sum(jfunsRing(p).*R4s,4);

                        cm = squeeze(sum(sum(jfun([pfit(1:13) 0]).*E4./C,2),1));
                        cmCore = squeeze(sum(sum(jfun([pfit(1) 0 0 pfit(4:13) 0]).*E4./C,2),1));
                        cmBeam = squeeze(sum(sum(jfunBeam(pfit).*E4./C,2),1));
                        cmA = squeeze(sum(sum(jfunA(pfit).*E4./C,2),1));
                        cmRing = squeeze(sum(sum(jfunRing(pfit).*E4./C,2),1));
                        cmBG = squeeze(sum(sum(repmat(pfit(15),size(E4)),2),1));

                        jqm4 = jfun(pfit);
                        Jqm = squeeze(sum(sum(jfun(pfit).*Omega4,2),1));
                        JqmCore = squeeze(sum(sum(jfun([pfit(1) 0 0 pfit(4:13) 0]).*Omega4,2),1));
                        JqmA    = squeeze(sum(sum(jfunA(pfit).*Omega4,2),1));
                        JqmBeam = squeeze(sum(sum(jfunBeam(pfit).*Omega4,2),1));
                        JqmRing = squeeze(sum(sum(jfunRing(pfit).*Omega4,2),1));
                        
                        JqmBG = squeeze(sum(sum(pfit(15)*C./E4.*Omega4,2),1));
                        Jqp = squeeze(sum(sum(jq.*Omega4,2),1));
                        sigJqp = sqrt(squeeze(sum(sum(sigjq.^2.*Omega4.^2,2),1)));

                        pp = pfit; 
                        jEp = jq.*E4; jEp = jEp(:); S1 = log10(jEp)+1; S1(S1 <= 0) = 0.1;
                        jEf = jqm4.*E4;
                        jEfv = jEf(:); S2 = log10(jEfv)+1; S2(S2 <= 7.1) = 0.1;
                        
                        [thvswnSWIA,phvswnSWIA,~] = cart2sph(SWIAf.vswnSWIA(1,indf),SWIAf.vswnSWIA(2,indf),SWIAf.vswnSWIA(3,indf));
                        if thvswnSWIA < 0; thvswnSWIA = thvswnSWIA+2*pi; end

                        %Spherical Natural distributions
                        dazel = 0.5;
                        Enp1 = 10.^(2:0.01:4.3);
                        azp1 = 140:dazel:220;
                        elp1 = -45:dazel:45;
                        [azp,elp,Eqp] = ndgrid(azp1,elp1,Enp1);
                        [~,~,dEqp] = ndgrid(azp1,elp1,gradient(Enp1));
                        [vHnx,vHny,vHnz] = sph2cart(deg2rad(azp),deg2rad(elp),sqrt(2*Eqp*eV/mH));
                        vHnAbs = sqrt(vHnx.^2 + vHny.^2 + vHnz.^2);
                        
                        [vAnx,vAny,vAnz] = sph2cart(deg2rad(azp),deg2rad(elp),sqrt(2*2*Eqp*eV/mA));
                        vAnAbs = sqrt(vAnx.^2 + vAny.^2 + vAnz.^2);

                        clearvars funnvHpara funnvHperp funnvApara funnvAperp;
                        funnvHpara = @(p) ((vHnx-p(4)*1000).*Bxn + (vHny-p(5)*1000).*Byn + (vHnz-p(6)*1000).*Bzn);
                        funnvHperp = @(p) sqrt(((vHnx-p(4)*1000)-funnvHpara(p).*Bxn).^2 + ((vHny-p(5)*1000)-funnvHpara(p).*Byn).^2 + ((vHnz-p(6)*1000)-funnvHpara(p).*Bzn).^2);
                        
                        funnvApara = @(p) ((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn).*Bxn + (vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn).*Byn + (vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn).*Bzn);
                        funnvAperp = @(p) sqrt(((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn)-funnvApara(p).*Bxn).^2 + ((vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn)-funnvApara(p).*Byn).^2 + ((vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn)-funnvApara(p).*Bzn).^2);
                        
                        funnvHparaRing = @(p) ((vHnx-Vrx*1000).*Bxn + (vHny-Vry*1000).*Byn + (vHnz-Vrz*1000).*Bzn);
                        funnvHperpRing = @(p) sqrt(((vHnx-Vrx*1000)-funnvHparaRing(p).*Bxn).^2 + ((vHny-Vry*1000)-funnvHparaRing(p).*Byn).^2 + ((vHnz-Vrz*1000)-funnvHparaRing(p).*Bzn).^2);
                    
                        jfunnCore = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(1).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(10))/(3*p(10)^(2/3)))^(3/2)*sqrt(2*k*p(7)/mH)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(10))/(3*sqrt(2*k*p(7)/mH)^2)*(funnvHpara(p).^2 + 1/p(10)*funnvHperp(p).^2)).^(-p(12)-1);
                        jfunnBeam = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(3).*(pi*sqrt(2*k*p(7)/mH)^2).^(-3/2).*(p(12)-3/2).^(-3/2).*gamma(p(12)+1)./gamma(p(12)-1/2).*(1 + ((vHnx-(p(4)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnxSWIA)).^2 + (vHny-(p(5)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnySWIA)).^2 + (vHnz-(p(6)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnzSWIA)).^2)./((p(12)-3/2).*sqrt(2*k*p(7)/mH).^2)).^(-p(12)-1);
                        jfunnA    = @(p) 2*1E-4*vAnAbs.^2./mA.*eV.*p(2).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(11))/(3*p(11)^(2/3)))^(3/2)*sqrt(2*k*p(9)*p(7)/mA)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(11))/(3*sqrt(2*k*p(9)*p(7)/mA)^2)*(funnvApara(p).^2 + 1/p(11)*funnvAperp(p).^2)).^(-p(12)-1);
                        jfunnRing = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(14)*(pi^(3/2)*sqrt(2*k*Tring/mH)^3*funA(p))^(-1).*exp(-(funnvHparaRing(p).^2)/(sqrt(2*k*Tring/mH))^2).*exp(-(funnvHperpRing(p)-(VrAbs*1000)).^2/(sqrt(2*k*Tring/mH))^2);
                        
                        jqns = jfunnCore(pfit) + jfunnBeam(pfit) + jfunnA(pfit) + jfunnRing(pfit);
                        jEns = jqns.*Eqp; 
                        
                        Omeganp = deg2rad(dazel).^2.*cosd(elp);
                        jEnsCoreOmni = squeeze(sum(sum(jfunnCore(pfit).*Eqp.*Omeganp,2),1)); 
                        jEnsBeamOmni = squeeze(sum(sum(jfunnBeam(pfit).*Eqp.*Omeganp,2),1));
                        jEnsAOmni = squeeze(sum(sum(jfunnA(pfit).*Eqp.*Omeganp,2),1));
                        jEnsRingOmni = squeeze(sum(sum(jfunnRing(pfit).*Eqp.*Omeganp,2),1));                        

                        %Cartesian Natural distributions
                        dvxp = -500; %Shift range in the vx-direction
                        dvzp = +0; %Shift range in the vz-direction
                        fEp = 1;%Multiplier for 1D-plots' energy range
                        dvxyz = 5; %Resolution in velocity-space [km/s]
                        [vHnx,vHny,vHnz] = meshgrid((-1000:dvxyz:-200)+dvxp,-400:dvxyz:400,(-400:dvxyz:400)+dvzp);
                        vHnx = 1000*vHnx; vHny = 1000*vHny; vHnz = 1000*vHnz;

                        vHnAbs = sqrt(vHnx.^2 + vHny.^2 + vHnz.^2);
                        En4 = mH*vHnAbs.^2/2/eV;

                        vAnx = vHnx/sqrt(2); vAny = vHny/sqrt(2); vAnz = vHnz/sqrt(2);
                        vAnAbs = sqrt(vAnx.^2 + vAny.^2 + vAnz.^2);

                        clearvars funnvHpara funnvHperp funnvApara funnvAperp funnvHparaRing funnvHperpRing jfunnCore jfunnBeam jfunnA jfunnRing;
                        funnvHpara = @(p) ((vHnx-p(4)*1000).*Bxn + (vHny-p(5)*1000).*Byn + (vHnz-p(6)*1000).*Bzn);
                        funnvHperp = @(p) sqrt(((vHnx-p(4)*1000)-funnvHpara(p).*Bxn).^2 + ((vHny-p(5)*1000)-funnvHpara(p).*Byn).^2 + ((vHnz-p(6)*1000)-funnvHpara(p).*Bzn).^2);
                        
                        funnvApara = @(p) ((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn).*Bxn + (vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn).*Byn + (vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn).*Bzn);
                        funnvAperp = @(p) sqrt(((vAnx-p(4)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bxn)-funnvApara(p).*Bxn).^2 + ((vAny-p(5)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Byn)-funnvApara(p).*Byn).^2 + ((vAnz-p(6)*1000-p(8)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*Bzn)-funnvApara(p).*Bzn).^2);
                        
                        funnvHparaRing = @(p) ((vHnx-Vrx*1000).*Bxn + (vHny-Vry*1000).*Byn + (vHnz-Vrz*1000).*Bzn);
                        funnvHperpRing = @(p) sqrt(((vHnx-Vrx*1000)-funnvHparaRing(p).*Bxn).^2 + ((vHny-Vry*1000)-funnvHparaRing(p).*Byn).^2 + ((vHnz-Vrz*1000)-funnvHparaRing(p).*Bzn).^2);
                    
                        jfunnCore = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(1).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(10))/(3*p(10)^(2/3)))^(3/2)*sqrt(2*k*p(7)/mH)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(10))/(3*sqrt(2*k*p(7)/mH)^2)*(funnvHpara(p).^2 + 1/p(10)*funnvHperp(p).^2)).^(-p(12)-1);
                        jfunnBeam = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(3).*(pi*sqrt(2*k*p(7)/mH)^2).^(-3/2).*(p(12)-3/2).^(-3/2).*gamma(p(12)+1)./gamma(p(12)-1/2).*(1 + ((vHnx-(p(4)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnxSWIA)).^2 + (vHny-(p(5)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnySWIA)).^2 + (vHnz-(p(6)*1000+p(13)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2)))*vnzSWIA)).^2)./((p(12)-3/2).*sqrt(2*k*p(7)/mH).^2)).^(-p(12)-1);
                        jfunnA    = @(p) 2*1E-4*vAnAbs.^2./mA.*eV.*p(2).*pi^(-3/2)*(p(12)-3/2)^(-3/2)*gamma(p(12)+1)/gamma(p(12)-0.5)*((1+2*p(11))/(3*p(11)^(2/3)))^(3/2)*sqrt(2*k*p(9)*p(7)/mA)^(-3).*(1 + 1/(p(12)-3/2)*(1+2*p(11))/(3*sqrt(2*k*p(9)*p(7)/mA)^2)*(funnvApara(p).^2 + 1/p(11)*funnvAperp(p).^2)).^(-p(12)-1);
                        jfunnRing = @(p)   1E-4*vHnAbs.^2./mH.*eV.*p(14)*(pi^(3/2)*sqrt(2*k*Tring/mH)^3*funA(p))^(-1).*exp(-(funnvHparaRing(p).^2)/(sqrt(2*k*Tring/mH))^2).*exp(-(funnvHperpRing(p)-(VrAbs*1000)).^2/(sqrt(2*k*Tring/mH))^2);
                        
                        jqn4 = jfunnCore(pfit) + jfunnBeam(pfit) + jfunnA(pfit) + jfunnRing(pfit);
                        jEnf = jqn4.*En4; 
                        
                    end
            end
            % Fit the windowed coarse distribution in 1D if no fine distributions have been included
            if isempty(indsFine) || sum(indsFineIncl) == 0 

                tic;
               
                clearvars p0 lb ub;
                p0(1) = 1.0*SWIAc.N(ii)*1E+6; %NH [m^-3]
                p0(2) = 0.03*SWIAc.N(ii)*1E+6;%NA [m^-3]
                p0(3) = 0.3*p0(1);
                p0(4) = -SWIAc.Vabs(ii);%Vabs
                p0(5) = 0.2*SWIAc.T(ii);%TH
                p0(6) = 1;%a |VA|/|VH|
                p0(7) = 4;%b TA/TH
                p0(8) = 0.01;
        
                lb(1) = 0.5*SWIAc.N(ii)*1E+6;
                lb(2) = 0;
                lb(3) = 0;
                lb(4) = -SWIAc.Vabs(ii)-50;
                lb(5) = 0.01*SWIAc.T(ii);
                lb(6) = 0.99;
                lb(7) = 0.8;
                lb(8) = 0;
        
                ub(1) = inf;
                ub(2) = 0.2*SWIAc.N(ii)*1E+6;
                ub(3) = 0.5*ub(1);
                ub(4) = -SWIAc.Vabs(ii)+50;
                ub(5) = 2*SWIAc.T(ii);
                ub(6) = 1.05;
                ub(7) = 10;
                ub(8) = 1000;

                %Skip empty distributions
                if sum(isnan([p0,lb,ub])) > 0
                    continue;
                end
                
                %Magnetic field vector
                Bxn = double(SWIAc.BxnSWIA(ii)); 
                Byn = double(SWIAc.BynSWIA(ii)); 
                Bzn = double(SWIAc.BznSWIA(ii));
                Babs = double(sqrt(SWIAc.BxSWIA(ii)^2 + SWIAc.BySWIA(ii)^2 + SWIAc.BzSWIA(ii)^2));
                
                %If B-field measurements not available, assume Parker spiral conditions
                if isnan(Babs)
                    Bxn = cosd(56);
                    Byn = sind(56);
                    Bzn = 0;
                    Babs = 2;
                end
        
                E4 = double(SWIAc.E4(:,:,:,ii)); 
                Omega4 = double(SWIAc.Omega4(:,:,:,ii));
                jq = double(SWIAc.j(:,:,:,ii));
                c = double(SWIAc.c(:,:,:,ii));
                C = double(SWIAc.C(:,:,:,ii));
                sigjq = double(C.*sqrt(c+1)./E4);

                c1 = squeeze(sum(sum(c,2),1)); E1 = SWIAc.E;
                Jq = squeeze(sum(sum(jq.*Omega4.*SWIAc.isWindow(:,:,:,ii),2),1));
                sigJq = sqrt(squeeze(sum(sum(sigjq.^2.*Omega4.^2.*SWIAc.isWindow(:,:,:,ii),2),1)));
        
                indSelE = squeeze(sum(sum(SWIAc.isWindow(:,:,:,ii),2),1)) > 1;
        
                Jq = Jq(indSelE);
                sigJq = sigJq(indSelE);
                E1 = E1(indSelE);
                c1 = c1(indSelE);

                %Cylindrically symmetric samples in velocity space
                dth = deg2rad(1);
        
                [thg,vHabs] = ndgrid(deg2rad(0:1:45)+pi,double(sqrt(2*E1*eV/mH)));
                [vHx,vHr,~] = sph2cart(thg,0,vHabs);
        
                [thg,vAabs] = ndgrid(deg2rad(0:1:45)+pi,double(sqrt(2*2*E1*eV/mA)));
                [vAx,vAr,~] = sph2cart(thg,0,vAabs);
        
                [~,E2] = ndgrid(deg2rad(0:1:45)+pi,double(E1));
                [~,C2] = ndgrid(deg2rad(0:1:45)+pi,double(squeeze(sum(sum(C(:,:,indSelE)/(4*pi),2),1))));

                %Fitting routine
                clearvars p pfit1dSlow pfit1dFast;
                
                JfunSlowBeam = @(p) squeeze(sum(2*pi*(thg-pi)*dth.*(1E-4.*vHabs.^2./mH.*eV.*p(1).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-p(4)*1000).^2 + vHr.^2)./(2*k*p(5)))... 
                                                             + 1E-4.*vHabs.^2./mH.*eV.*p(3).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-(p(4)*1000+abs(Bxn)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2))))).^2 + vHr.^2)./(2*k*p(5)))...
                                                           + 2*1E-4.*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(7)*p(5)))^(3/2).*exp(-mA.*((vAx-p(6)*p(4)*1000).^2 + vAr.^2)./(2*k*p(7)*p(5)))),1));%    + 2*1E-4*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(8)*p(6))).^(3/2).*exp(-mA.*((vAx-p(7).*p(3)*1000).^2 + (vAy-p(7).*p(4)*1000).^2 + (vAz-p(7).*p(5)*1000).^2)./(2*k*p(8)*p(6)));
                
                JfunFastBeam = @(p) squeeze(sum(2*pi*(thg-pi)*dth.*(1E-4.*vHabs.^2./mH.*eV.*p(1).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-p(4)*1000).^2 + vHr.^2)./(2*k*p(5)))... 
                                                             + 1E-4.*vHabs.^2./mH.*eV.*p(3).*(mH./(2*pi*k*p(5)))^(3/2).*exp(-mH.*((vHx-(p(4)*1000-abs(Bxn)*Babs*1E-9/sqrt(mu0*(mH*(p(1)+p(3)) + mA*p(2))))).^2 + vHr.^2)./(2*k*p(5)))...
                                                           + 2*1E-4.*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(7)*p(5)))^(3/2).*exp(-mA.*((vAx-p(6)*p(4)*1000).^2 + vAr.^2)./(2*k*p(7)*p(5)))),1));%    + 2*1E-4*vAabs.^2./mA.*eV.*p(2).*(mA./(2*pi*k*p(8)*p(6))).^(3/2).*exp(-mA.*((vAx-p(7).*p(3)*1000).^2 + (vAy-p(7).*p(4)*1000).^2 + (vAz-p(7).*p(5)*1000).^2)./(2*k*p(8)*p(6)));
                
                costfunSlowBeam = @(p) sum((Jq - JfunSlowBeam(p)').^2./sigJq.^2); %Chi-square cost function
                costfunFastBeam = @(p) sum((Jq - JfunFastBeam(p)').^2./sigJq.^2); %Chi-square cost function
                
                opts = optimoptions('fmincon','Display','off');
                [pfit1dSlow,chi2Slow,fitExitFlagSlow,~,~,~,Hslow] = fmincon(costfunSlowBeam,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
                [pfit1dFast,chi2Fast,fitExitFlagFast,~,~,~,Hfast] = fmincon(costfunFastBeam,double(p0),[],[],[],[],double(lb),double(ub),[],opts); %Bound GN-ish optimization
                
                chi2rSlow = chi2Slow/(length(Jq)-length(pfit1dSlow));
                chi2rFast = chi2Fast/(length(Jq)-length(pfit1dFast));
                
                if chi2rFast < chi2rSlow
                    pfit1d = pfit1dFast;
                    chi2r = chi2rFast;
                    fitExitFlag = fitExitFlagFast;
                    H = Hfast;
                    Jfun = JfunFastBeam;
                    isFastBeam = true;
                else
                    pfit1d = pfit1dSlow;
                    chi2r = chi2rSlow;
                    fitExitFlag = fitExitFlagSlow;
                    H = Hslow;
                    Jfun = JfunSlowBeam;
                    isFastBeam = false;
                end
                
                sigpfit1d = sqrt(diag(inv(H/2)));
        
                if sum(imag(sigpfit1d) > 0)
                    warning('Bad FIM from 1D fit, solution likely at boundary');
                    clearvars sigpfit1d;
                    sigpfit1d = inf(1,8);
                end
                
                VHabs = (pfit1d(1)*pfit1d(4) + pfit1d(3)*(pfit1d(4) + 1E-3*abs(Bxn)*Babs*1E-9/sqrt(mu0*(mH*(pfit1d(1)+pfit1d(3)) + mA*pfit1d(2)))))/(pfit1d(1) + pfit1d(3));
                VAabs =  pfit1d(6)*pfit1d(4);
                VHabs = abs(VHabs);
                VAabs = abs(VAabs);
                TH = pfit1d(5);
                TA = pfit1d(7)*pfit1d(5);

                sigNH = sqrt(sigpfit1d(1)^2+sigpfit1d(3)^2)*1E-6;
                sigNA = sigpfit1d(2)*1E-6;
                sigTH = sigpfit1d(5);
                sigTA = sqrt(TA^2*((sigpfit1d(7)/pfit1d(7))^2 + (sigpfit1d(5)/pfit1d(5))^2));
                
                Jq = single(zeros(1,48));
                Jq(indSelE) = Jfun(pfit1d);
                
                pfit = NaN(1,15); pfit(1:length(pfit1d)) = pfit1d;

                %Analyze whether the distribution is likely undisturbed solar wind
                isSW = true;
                VthH = sqrt(2*k*TH/mH); 
                
                %Thermal speed too high (likely sheath)
                if VthH/(1000*VHabs) > 0.3
                    isSW = false;
                end
                
                %Concatenate common structure
                SWIA.time = cat(1,SWIA.time,SWIAc.time(ii));
                SWIA.NH = cat(1,SWIA.NH,(pfit1d(1)+pfit1d(3))*1E-6);
                SWIA.NA = cat(1,SWIA.NA,pfit1d(2)*1E-6);
                SWIA.NHcore    = cat(1,SWIA.NHcore,NaN);
                SWIA.NHbeam    = cat(1,SWIA.NHbeam,NaN);
                SWIA.VHx       = cat(1,SWIA.VHx,NaN);
                SWIA.VHy       = cat(1,SWIA.VHy,NaN);
                SWIA.VHz       = cat(1,SWIA.VHz,NaN);
                SWIA.VHabs     = cat(1,SWIA.VHabs,VHabs);
                SWIA.VAx       = cat(1,SWIA.VAx,NaN);
                SWIA.VAy       = cat(1,SWIA.VAy,NaN);
                SWIA.VAz       = cat(1,SWIA.VAz,NaN);
                SWIA.VAabs     = cat(1,SWIA.VAabs,VAabs);
                SWIA.THcore    = cat(1,SWIA.THcore,TH);
                SWIA.TA        = cat(1,SWIA.TA,TA);
                SWIA.TH     = cat(1,SWIA.TH,NaN);
                SWIA.THperp = cat(1,SWIA.THperp,NaN);
                SWIA.THpara = cat(1,SWIA.THpara,NaN);
                SWIA.alphaHcore = cat(1,SWIA.alphaHcore,NaN);
                SWIA.alphaA     = cat(1,SWIA.alphaA,NaN);
                SWIA.kappa     = cat(1,SWIA.kappa,NaN);

                SWIA.VHcorex = cat(1,SWIA.VHcorex,NaN);
                SWIA.VHcorey = cat(1,SWIA.VHcorey,NaN);
                SWIA.VHcorez = cat(1,SWIA.VHcorez,NaN);
                SWIA.dVHbeamx = cat(1,SWIA.dVHbeamx,NaN);
                SWIA.dVHbeamy = cat(1,SWIA.dVHbeamy,NaN);
                SWIA.dVHbeamz = cat(1,SWIA.dVHbeamz,NaN);
                
                SWIA.Bx        = cat(1,SWIA.Bx,SWIAc.Bx(ii));
                SWIA.By        = cat(1,SWIA.By,SWIAc.By(ii));
                SWIA.Bz        = cat(1,SWIA.Bz,SWIAc.Bz(ii));

                SWIA.BxSWIA    = cat(1,SWIA.BxSWIA,SWIAc.BxSWIA(ii));
                SWIA.BySWIA    = cat(1,SWIA.BySWIA,SWIAc.BySWIA(ii));
                SWIA.BzSWIA    = cat(1,SWIA.BzSWIA,SWIAc.BzSWIA(ii));

                SWIA.vnxSWIA   = cat(1,SWIA.vnxSWIA,NaN);
                SWIA.vnySWIA   = cat(1,SWIA.vnySWIA,NaN);
                SWIA.vnzSWIA   = cat(1,SWIA.vnzSWIA,NaN);

                SWIA.Valf   = cat(1,SWIA.Valf,NaN);

                SWIA.En           = cat(1,SWIA.En,10.^linspace(1,4,numelsE));
                SWIA.jEnCoreOmni = cat(1,SWIA.jEnCoreOmni,NaN(1,numelsE));
                SWIA.jEnBeamOmni = cat(1,SWIA.jEnBeamOmni,NaN(1,numelsE));
                SWIA.jEnAomni    = cat(1,SWIA.jEnAomni,NaN(1,numelsE));
                SWIA.jEnRingOmni = cat(1,SWIA.jEnRingOmni,NaN(1,numelsE));

                SWIA.sigNH     = cat(1,SWIA.sigNH,sigNH);
                SWIA.sigNA     = cat(1,SWIA.sigNA,sigNA);

                SWIA.sigNHcore     = cat(1,SWIA.sigNHcore,NaN);
                SWIA.sigNHbeam     = cat(1,SWIA.sigNHbeam,NaN);

                SWIA.sigVHx       = cat(1,SWIA.sigVHx,NaN);
                SWIA.sigVHy       = cat(1,SWIA.sigVHy,NaN);
                SWIA.sigVHz       = cat(1,SWIA.sigVHz,NaN);
                SWIA.sigVAx       = cat(1,SWIA.sigVAx,NaN);
                SWIA.sigVAy       = cat(1,SWIA.sigVAy,NaN);
                SWIA.sigVAz       = cat(1,SWIA.sigVAz,NaN);

                SWIA.sigVHcorex    = cat(1,SWIA.sigVHcorex,NaN);
                SWIA.sigVHcorey    = cat(1,SWIA.sigVHcorey,NaN);
                SWIA.sigVHcorez    = cat(1,SWIA.sigVHcorez,NaN);
                SWIA.sigdVHbeamx   = cat(1,SWIA.sigdVHbeamx,NaN);
                SWIA.sigdVHbeamy   = cat(1,SWIA.sigdVHbeamy,NaN);
                SWIA.sigdVHbeamz   = cat(1,SWIA.sigdVHbeamz,NaN);
                SWIA.sigTHcore     = cat(1,SWIA.sigTHcore,sigTH);
                SWIA.sigTA         = cat(1,SWIA.sigTA,sigTA);
                SWIA.sigAlphaHcore = cat(1,SWIA.sigAlphaHcore,NaN);
                SWIA.sigAlphaA     = cat(1,SWIA.sigAlphaA,NaN);
                SWIA.sigKappa      = cat(1,SWIA.sigKappa,NaN);

                SWIA.Eq        = cat(1,SWIA.Eq,SWIAc.E');
                SWIA.Jq        = cat(1,SWIA.Jq,Jq);
                SWIA.pfit      = cat(1,SWIA.pfit,pfit);
                SWIA.frac95    = cat(1,SWIA.frac95,NaN);
                SWIA.chi2rApprox = cat(1,SWIA.chi2rApprox,NaN);

                SWIA.indrSel   = cat(1,SWIA.indrSel,NaN);
                SWIA.atten_state = cat(1,SWIA.atten_state,SWIAc.atten_state(ii));
                
                SWIA.Mswia2mso = cat(1,SWIA.Mswia2mso,shiftdim(SWIAc.Mswia2msoAll(:,:,ii),-1));
                SWIA.Mmso2swia = cat(1,SWIA.Mmso2swia,shiftdim(SWIAc.Mmso2swiaAll(:,:,ii),-1));
                SWIA.isSW      = cat(1,SWIA.isSW,isSW);
                SWIA.isDraped  = cat(1,SWIA.isDraped,SWIAc.isDraped(ii));
                SWIA.isIMF     = cat(1,SWIA.isIMF,SWIAc.isIMF(ii));
                SWIA.isIono    = cat(1,SWIA.isIono,SWIAc.isIono(ii));
                SWIA.isFine    = cat(1,SWIA.isFine,false);
                SWIA.isFS      = cat(1,SWIA.isFS,SWIAc.isFS(ii));
                SWIA.isSvy     = cat(1,SWIA.isSvy,SWIAc.grouping(ii) == 1);
                SWIA.fitExitFlag = cat(1,SWIA.fitExitFlag,fitExitFlag);
                
                SWIA.xSC     = cat(1,SWIA.xSC,SWIAc.xSC(ii));
                SWIA.ySC     = cat(1,SWIA.ySC,SWIAc.ySC(ii));
                SWIA.zSC     = cat(1,SWIA.zSC,SWIAc.zSC(ii));

                SWIA.dtRunTime = cat(1,SWIA.dtRunTime,toc);

                indsFineLog = cat(1,indsFineLog,NaN);
            end
        else %If not identified as solar wind distribution

            %Concatenate common structure
            tic;
            if ~isnan(SWIAc.NHfit(ii))
                SWIA.time      = cat(1,SWIA.time,SWIAc.time(ii));
                SWIA.NH        = cat(1,SWIA.NH,SWIAc.NHfit(ii));
                SWIA.NA        = cat(1,SWIA.NA,SWIAc.NAfit(ii));
                SWIA.NHcore    = cat(1,SWIA.NHcore,NaN);
                SWIA.NHbeam    = cat(1,SWIA.NHbeam,NaN);
                SWIA.VHx       = cat(1,SWIA.VHx,NaN);
                SWIA.VHy       = cat(1,SWIA.VHy,NaN);
                SWIA.VHz       = cat(1,SWIA.VHz,NaN);
                SWIA.VHabs     = cat(1,SWIA.VHabs,SWIAc.VHabsFit(ii));
                SWIA.VAx       = cat(1,SWIA.VAx,NaN);
                SWIA.VAy       = cat(1,SWIA.VAy,NaN);
                SWIA.VAz       = cat(1,SWIA.VAz,NaN);
                SWIA.VAabs     = cat(1,SWIA.VAabs,SWIAc.VAabsFit(ii));
                SWIA.TH        = cat(1,SWIA.TH,SWIAc.THfit(ii));
                SWIA.TA        = cat(1,SWIA.TA,SWIAc.TAfit(ii));
                SWIA.THcore     = cat(1,SWIA.THcore,NaN);
                SWIA.THperp = cat(1,SWIA.THperp,NaN);
                SWIA.THpara = cat(1,SWIA.THpara,NaN);
                SWIA.alphaHcore = cat(1,SWIA.alphaHcore,NaN);
                SWIA.alphaA     = cat(1,SWIA.alphaA,NaN);
                SWIA.kappa     = cat(1,SWIA.kappa,NaN);

                SWIA.VHcorex = cat(1,SWIA.VHcorex,NaN);
                SWIA.VHcorey = cat(1,SWIA.VHcorey,NaN);
                SWIA.VHcorez = cat(1,SWIA.VHcorez,NaN);
                SWIA.dVHbeamx = cat(1,SWIA.dVHbeamx,NaN);
                SWIA.dVHbeamy = cat(1,SWIA.dVHbeamy,NaN);
                SWIA.dVHbeamz = cat(1,SWIA.dVHbeamz,NaN);
                
                SWIA.Bx        = cat(1,SWIA.Bx,SWIAc.Bx(ii));
                SWIA.By        = cat(1,SWIA.By,SWIAc.By(ii));
                SWIA.Bz        = cat(1,SWIA.Bz,SWIAc.Bz(ii));
    
                SWIA.BxSWIA    = cat(1,SWIA.BxSWIA,SWIAc.BxSWIA(ii));
                SWIA.BySWIA    = cat(1,SWIA.BySWIA,SWIAc.BySWIA(ii));
                SWIA.BzSWIA    = cat(1,SWIA.BzSWIA,SWIAc.BzSWIA(ii));
    
                SWIA.vnxSWIA   = cat(1,SWIA.vnxSWIA,NaN);
                SWIA.vnySWIA   = cat(1,SWIA.vnySWIA,NaN);
                SWIA.vnzSWIA   = cat(1,SWIA.vnzSWIA,NaN);

                SWIA.Valf   = cat(1,SWIA.Valf,NaN);

                SWIA.En           = cat(1,SWIA.En,10.^linspace(1,4,numelsE));
                SWIA.jEnCoreOmni = cat(1,SWIA.jEnCoreOmni,NaN(1,numelsE));
                SWIA.jEnBeamOmni = cat(1,SWIA.jEnBeamOmni,NaN(1,numelsE));
                SWIA.jEnAomni    = cat(1,SWIA.jEnAomni,NaN(1,numelsE));
                SWIA.jEnRingOmni = cat(1,SWIA.jEnRingOmni,NaN(1,numelsE));
    
                SWIA.sigNH     = cat(1,SWIA.sigNH,NaN);
                SWIA.sigNA     = cat(1,SWIA.sigNA,NaN);

                SWIA.sigNHcore     = cat(1,SWIA.sigNHcore,NaN);
                SWIA.sigNHbeam     = cat(1,SWIA.sigNHbeam,NaN);

                SWIA.sigVHx     = cat(1,SWIA.sigVHx,NaN);
                SWIA.sigVHy     = cat(1,SWIA.sigVHy,NaN);
                SWIA.sigVHz     = cat(1,SWIA.sigVHz,NaN);
                SWIA.sigVAx     = cat(1,SWIA.sigVAx,NaN);
                SWIA.sigVAy     = cat(1,SWIA.sigVAy,NaN);
                SWIA.sigVAz     = cat(1,SWIA.sigVAz,NaN);

                SWIA.sigVHcorex = cat(1,SWIA.sigVHcorex,NaN);
                SWIA.sigVHcorey = cat(1,SWIA.sigVHcorey,NaN);
                SWIA.sigVHcorez = cat(1,SWIA.sigVHcorez,NaN);
                SWIA.sigdVHbeamx = cat(1,SWIA.sigdVHbeamx,NaN);
                SWIA.sigdVHbeamy = cat(1,SWIA.sigdVHbeamy,NaN);
                SWIA.sigdVHbeamz = cat(1,SWIA.sigdVHbeamz,NaN);

                SWIA.sigTA     = cat(1,SWIA.sigTA,NaN);
                SWIA.sigTHcore = cat(1,SWIA.sigTHcore,NaN);
                SWIA.sigAlphaHcore  = cat(1,SWIA.sigAlphaHcore,NaN);
                SWIA.sigAlphaA  = cat(1,SWIA.sigAlphaA,NaN);
                SWIA.sigKappa  = cat(1,SWIA.sigKappa,NaN);
    
                SWIA.Eq        = cat(1,SWIA.Eq,SWIAc.E');
                SWIA.Jq        = cat(1,SWIA.Jq,SWIAc.JqFit(:,ii)');
                SWIA.pfit      = cat(1,SWIA.pfit,NaN(1,15));
                SWIA.frac95    = cat(1,SWIA.frac95,NaN);
                SWIA.chi2rApprox = cat(1,SWIA.chi2rApprox,NaN);

                SWIA.indrSel   = cat(1,SWIA.indrSel,NaN);
                SWIA.atten_state = cat(1,SWIA.atten_state,SWIAc.atten_state(ii));
                
                SWIA.Mswia2mso = cat(1, SWIA.Mswia2mso,shiftdim(SWIAc.Mswia2msoAll(:,:,ii),-1));
                SWIA.Mmso2swia = cat(1, SWIA.Mmso2swia,shiftdim(SWIAc.Mmso2swiaAll(:,:,ii),-1));
                SWIA.isSW      = cat(1,SWIA.isSW,SWIAc.isSW(ii));
                SWIA.isDraped  = cat(1,SWIA.isDraped,SWIAc.isDraped(ii));
                SWIA.isIMF     = cat(1,SWIA.isIMF,SWIAc.isIMF(ii));
                SWIA.isIono    = cat(1,SWIA.isIono,SWIAc.isIono(ii));
                SWIA.isFine    = cat(1,SWIA.isFine,false);
                SWIA.isFS      = cat(1,SWIA.isFS,SWIAc.isFS(ii));
                SWIA.isSvy     = cat(1,SWIA.isSvy,SWIAc.grouping(ii) == 1);
                SWIA.fitExitFlag = cat(1,SWIA.fitExitFlag,NaN);
                
                SWIA.xSC     = cat(1,SWIA.xSC,SWIAc.xSC(ii));
                SWIA.ySC     = cat(1,SWIA.ySC,SWIAc.ySC(ii));
                SWIA.zSC     = cat(1,SWIA.zSC,SWIAc.zSC(ii));
    
                SWIA.dtRunTime = cat(1,SWIA.dtRunTime,toc);

                indsFineLog = cat(1,indsFineLog,NaN);
            end
            
        end
    end

    SWIA.Babs = sqrt(SWIA.Bx.^2 + SWIA.By.^2 + SWIA.Bz.^2);
    
    %Calculate parallel/perpendicular temperatures
    SWIA.THcorePara = 3./(1+2*SWIA.alphaHcore).*SWIA.THcore;
    SWIA.THcorePerp = 3.*SWIA.alphaHcore./(1+2*SWIA.alphaHcore).*SWIA.THcore;

    SWIA.TApara = 3./(1+2*SWIA.alphaA).*SWIA.TA;
    SWIA.TAperp = 3.*SWIA.alphaA./(1+2*SWIA.alphaA).*SWIA.TA;

    SWIA.sigTHcorePara = sqrt((3./(1+2*SWIA.alphaHcore)).^2.*SWIA.sigTHcore.^2 + (6./(1+2*SWIA.alphaHcore).^2.*SWIA.THcore).^2.*SWIA.sigAlphaHcore.^2);
    SWIA.sigTHcorePerp = sqrt((3*SWIA.alphaHcore./(1+2*SWIA.alphaHcore)).^2.*SWIA.sigTHcore.^2 + (3.*SWIA.THcore./(1+2*SWIA.alphaHcore) - 6.*SWIA.THcore.*SWIA.alphaHcore./(1+2*SWIA.alphaHcore).^2).^2.*SWIA.sigAlphaHcore.^2);
    
    SWIA.sigTApara = sqrt((3./(1+2*SWIA.alphaA)).^2.*SWIA.sigTA.^2 + (6./(1+2*SWIA.alphaA).^2.*SWIA.TA).^2.*SWIA.sigAlphaA.^2);
    SWIA.sigTAperp = sqrt((3*SWIA.alphaA./(1+2*SWIA.alphaA)).^2.*SWIA.sigTA.^2 + (3.*SWIA.TA./(1+2*SWIA.alphaA) - 6.*SWIA.TA.*SWIA.alphaA./(1+2*SWIA.alphaA).^2).^2.*SWIA.sigAlphaA.^2);
    
    SWIA.alphaH = SWIA.THperp./SWIA.THpara;

    %Include the reponse functions in the saved file
    SWIA.SWIAro = SWIAro;
    SWIA.SWIArc = SWIArc;

    %Save the model functions (if used)
    try SWIA.funvHparaStr = char(funvHpara); catch; end
    try SWIA.funvHperpStr = char(funvHperp); catch; end
    try SWIA.funvAparaStr = char(funvApara); catch; end
    try SWIA.funvAperpStr = char(funvAperp); catch; end
    try SWIA.funAstr = char(funA); catch; end
    try SWIA.jfunsStr = char(jfuns); catch; end
    try SWIA.jfunStr = char(jfun); catch; end
    try SWIA.JfunStr = char(Jfun); catch; end
       
    %Check for errors in energy table
    V = SWIA.VHabs;
    %SWIAf.indBad = false(size(SWIAf.time));
    clearvars avgVSW stdVSW errLog
    for ii = 1:length(SWIA.time)
        %indw = find((abs(SWIA.time(ii)-SWIA.time) < 1/(60*24)).*((~SWIA.isSW + ~SWIA.isFine) > 0)); %Indices in moving window
        indw = find((abs(SWIA.time(ii)-SWIA.time) < 1/(60*24)).*SWIA.isSW.*SWIA.isFine); %Indices in moving window
        indw(indw == ii) = [];
        if length(indw) > 10 %Generally 15-30 distributions fit in a 2 min wide window
            
            avgV(ii) = mean(V(indw));
            stdV(ii) = std(V(indw));
            errLog(ii) = abs(V(ii)-avgV(ii))/stdV(ii);
            if abs(V(ii)-avgV(ii))/avgV(ii) > 0.1 && abs(V(ii)-avgV(ii))/stdV(ii) > 20
                disp(['   %% Found outlier in SWIA fine distribution #',num2str(ii),' on ',datestr(SWIA.time(ii),'yyyy-mm-dd HH:MM:SS'),' (',num2str(errLog(ii)),' sigma), flagging and removing from isSW']);
                
                SWIA.isSW(ii) = false; %Flag as a bad distribution
            end
        end
    end
    
    %Check for non-unique time stamps
    if length(unique(SWIA.time)) ~= length(SWIA.time)
        error('This orbit contains non-unique time-stamps, the checks didn''t work');
    end

    %Make sure moments cannot have mutually exclusive domains
    SWIA.isDraped(SWIA.isSW) = false;
    SWIA.isIono(SWIA.isSW) = false;

    if isempty(SWIA.time)
        warning(['Orbit #',num2str(OrbNum(indOrb)),' is empty, skipping save']);
        continue;
    end

    save([saveLoc,'mvn_swi_sw_',sprintf('%05d',OrbNum(indOrb)),'_v07_r01.mat'],'-struct','SWIA');
    
    disp(['  %% Completed MAVEN orbit #',num2str(OrbNum(indOrb)),':']);
    disp(['   %% Number of solar wind moments: ',num2str(sum(SWIA.isSW)),'/',num2str(length(SWIA.time)),' ']);
    disp(['   %% Response functions used: ',sprintf('%.0f ',sum(SWIA.indrSel == 0:length(SWIArc)))]);
    disp(['   %% Mean runtime/fine distribution: ',sprintf('%.1f',nanmean(SWIA.dtRunTime(SWIA.isFine))),' +- ',sprintf('%.1f',nanstd(SWIA.dtRunTime(SWIA.isFine))),' s']);

    %
    runTime = seconds(toc(tstartOrbit));
    runTime.Format = 'hh:mm:ss';
    disp(['  %% Total orbit runtime: ',char(runTime),' (h min s)']);

end
