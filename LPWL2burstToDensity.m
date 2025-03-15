%LPWL2burstToDensity

% Script to analyze MAVEN/LPW burst electric field data. Performs an FFT on
% the waveform and finds any prominent single-frequency peaks. The peak 
% frequency is taken as the plasma frequency and used to calculate the
% electron density.
% 
% Version history:
% 1.0 - Initial implementation
%
% Author: Robin Ramstad (2021-2025, LASP)

clear;

%Physical constants
 amu = 1.6605E-27; %[kg]
  me = 5.4859E-4*amu; %Electron mass
eps0 = 8.854187E-12; %Electric permittivity [F/m]
   e = 1.602E-19;% Elementary charge [C]

%Settings
doplot = 1;
makeCompleteFile = 0;
version = 'v1.0';
Orbits = 1:22318; %2014 to Oct 2024

%Load list and times of MAVEN orbits
load('MAVENorbits.mat'); 

%Load MAVEN SPK SPICE files
loadMAVENspiceV2spkOnly;

%% Go through orbit-by-orbit
indOrbs = find(ismember(OrbNum,Orbits));
for indNum = 1:length(indOrbs) 
    indOrb = indOrbs(indNum);
    
    filename = ['mvn_lpw_bursthf_orb',sprintf('%05d',OrbNum(indOrb)),'.mat'];
    
    disp(['%% Loading: ',filename]);
    try LPW = load([L2matLoc,filename]);
    catch cerror
        warning([' Could not load ',filename,' skipping - Cause: ',cerror.message]);
        continue;
    end
    
   %Identify burst groups
    dtt2000 = diff(double(LPW.tt2000));
    indBGroup = NaN(length(dtt2000),1);
    nBgroup = 1;
    for ii = 1:length(dtt2000)-1
        if dtt2000(ii) < 500
            indBGroup(ii) = nBgroup;
        else
            nBgroup = nBgroup + 1;
        end
    end
    
    %Initialize arrays
        time =        NaN(max(nBgroup),1);
          ne = single(NaN(max(nBgroup),1));
       signe = single(NaN(max(nBgroup),1));
           P = single(NaN(max(nBgroup),2^11));
       Epeak = single(NaN(max(nBgroup),1));
        prom = single(NaN(max(nBgroup),1));
       wpeak = single(NaN(max(nBgroup),1));
       fpeak = single(NaN(max(nBgroup),1));
    qualPeak = single(NaN(max(nBgroup),1));

    %FFT
    if doplot == 1; figure(1); clf; end
    dtInt = 500; %Interpolation time-step [ns]
    T = dtInt*1E-9;  
    Fs = 1/T; 
    n = 2^12; 
    f = 0:(Fs/n):(Fs/2-Fs/n); %Frequency domain (1-sided)
    for ii = 1:max(nBgroup)
        ind = indBGroup == ii;

        time(ii) = datenum(mean(LPW.timeDT(ind)));

        %Interpolate to remove FFT artifacts caused by gaps
        tt2000int = double(min(LPW.tt2000(ind))+dtInt:dtInt:max(LPW.tt2000(ind))-dtInt)';
        Eint = interp1(double(LPW.tt2000(ind)),LPW.E(ind),tt2000int);

        Eint = Eint-mean(Eint); %Remove non-zero offset

        L = length(tt2000int); 
        Y = fft(Eint,n);
        
        if sum(isnan(Y)) == length(Y)
            continue;
        end
        if sum(Y) == 0
            continue;
        end

        indSel = 1:n/2; %Select only 1-sided frequency domain
        P2 = abs(Y/n); %Double-sided power-spectrum
        P1 = P2(indSel); %Single-sided
        P1(2:end-1) = 2*P1(2:end-1); %Multiply by 2 to include the spectral power from the negative frequencies

        P(ii,:) = P1;

        %% Search for peaks
        indfsearch = ones(size(f)); 
        indfsearch((f > 8E+4).*(f < 1.1E+5) == 1) = 0;
        indfsearch((f > 2.5E+5).*(f < 3E+5) == 1) = 0;
        indfsearch(f < 3E+3) = 0;
        indfsearch = indfsearch == 1;
        fsearch = f(indfsearch);
        Psearch = P1(indfsearch);
        
        %Find the peak frequency
        [pks,fpks,w,p] = findpeaks(Psearch,fsearch);

        %Sort peak candidates by prominance
        [~,I] = sort(p,'descend');
        pks = pks(I);
        fpks = fpks(I);
        w = w(I);
        p = p(I);

        qualPeak(ii) = p(1)./p(2); %Peak quality index
        Epeak(ii) = pks(1);
        fpeak(ii) = fpks(1);
        prom(ii) = p(1);
        wpeak(ii) = w(1);
        
        omega = 2*pi*fpeak(ii);
        sigOmega = 2*pi*wpeak(ii);
        dnedomega = 2*me*eps0/e^2*omega;
        ne(ii) = me*eps0*omega^2/e^2*1E-6;
        signe(ii) = sqrt(dnedomega.^2.*sigOmega^2)*1E-6;
        
    end
    
    %% Derive locations using SPICE
    tstr = datestr(time);
    et = cspice_str2et(tstr); %Convert date string to ephemeris time
    [state, ltime] = cspice_spkezr('MAVEN', et, 'MAVEN_MSO', 'NONE', 'MARS');
    state = state';
    xMSO = single(state(:,1));
    yMSO = single(state(:,2));
    zMSO = single(state(:,3));
    
    %Get geographical coordinates
    [stateGEO,timeLocal] = cspice_spkezr('MAVEN', et, 'IAU_MARS', 'NONE', 'MARS');
    stateGEO = stateGEO';
    [~,altAreoid] = cspice_nearpt(state(:,1:3)',3396.9,3396.9,3376.097);
    [lon,lat,~] = cart2sph(stateGEO(:,1),stateGEO(:,2),stateGEO(:,3));
    lon = single(rad2deg(lon)); lat = single(rad2deg(lat));
    altAreoid = single(altAreoid)';
    
    %% Save the derived data in a file
    f = single(f);
    filename = ['mvn_lpw_hfFFTdens_orb',sprintf('%05d',OrbNum(indOrb)),'.mat'];
    disp([' %% Saving: ',filename]);
    save([L3matLoc,filename],'time','f','P','ne','signe','Epeak','prom','wpeak','fpeak','qualPeak','xMSO','yMSO','zMSO','lon','lat','altAreoid','version');
    
end
%% 
if makeCompleteFile == 1
   
    disp(['Making a complete file by concatenating all available files in ',L3matLoc]);
    
    list = dir(L3matLoc);
    
    %Initialize arrays
        time = double(NaN(0,1));
          ne = single(NaN(0,1));
       signe = single(NaN(0,1));
       PULF1 = single(NaN(0,1)); %<1 kHz
       PULF3 = single(NaN(0,1)); %<3 kHz
       PULF5 = single(NaN(0,1)); %<5 kHz
       Epeak = single(NaN(0,1));
        prom = single(NaN(0,1));
       wpeak = single(NaN(0,1));
       fpeak = single(NaN(0,1));
    qualPeak = single(NaN(0,1));
        xMSO = single(NaN(0,1));
        yMSO = single(NaN(0,1));
        zMSO = single(NaN(0,1));
         lon = single(NaN(0,1));
         lat = single(NaN(0,1));
   altAreoid = single(NaN(0,1));
         orb = single(NaN(0,1));
    for ii = 3:length(list)
        
        filename = list(ii).name;
        if ~strcmp(filename(1:17),'mvn_lpw_hfFFTdens')
            continue;
        end
        
        disp([' %% Concatenating: ',filename]);
        file = load([L3matLoc,filename]);
        
        df = gradient(file.f);
        
        if size(file.altAreoid,1) == 1
            file.altAreoid = file.altAreoid';
        end
        
         time = cat(1,time,file.time);
           ne = cat(1,ne,file.ne);
        signe = cat(1,signe,file.signe);
        PULF1 = cat(1,PULF1,sum(file.P(:,file.f < 1E+3).*repmat(df(file.f < 1E+3),[length(file.time) 1]),2));
        PULF3 = cat(1,PULF3,sum(file.P(:,file.f < 3E+3).*repmat(df(file.f < 3E+3),[length(file.time) 1]),2));
        PULF5 = cat(1,PULF5,sum(file.P(:,file.f < 5E+3).*repmat(df(file.f < 5E+3),[length(file.time) 1]),2));
        Epeak = cat(1,Epeak,file.Epeak);
        prom = cat(1,prom,file.Epeak);
        wpeak = cat(1,wpeak,file.wpeak);
        fpeak = cat(1,fpeak,file.fpeak);
        qualPeak = cat(1,qualPeak,file.qualPeak);
        xMSO = cat(1,xMSO,file.xMSO);
        yMSO = cat(1,yMSO,file.yMSO);
        zMSO = cat(1,zMSO,file.zMSO);
        lon = cat(1,lon,file.lon);
        lat = cat(1,lat,file.lat);
        altAreoid = cat(1,altAreoid,file.altAreoid);
        orb = cat(1,orb,single(str2double(filename(22:26)).*ones(size(file.time))));
    end
    
    f = file.f;

    disp(['Saving full concatenated file as ','LPWburstComplete.mat']);
    save([L3matLoc,'LPWburstComplete.mat'],'time','f','ne','signe','Epeak','prom','wpeak','fpeak','qualPeak','PULF1','PULF3','PULF5','xMSO','yMSO','zMSO','lon','lat','altAreoid','orb');
end