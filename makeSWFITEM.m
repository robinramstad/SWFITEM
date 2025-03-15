%makeSWFITEM
%
% Script to read files with fitted SWIA moments, reduce, concatenate, and 
% save as easily usable files.
%
% Author: Robin Ramstad (CU-LASP) 2025

clear;

disp(' % makeSWFITEM');

%Search paths
saveLoc = '.\SWFITEM\'; %Will save files in this folder

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
  k = 1.3806E-23;% [J/K] 8.6171E-5; %[eV/K] 
  q = 1.602176E-19; %Elementary charge [C] 
mu0 = 4*pi*1E-7; %Magnetic permeability in vacuum [H/m]

%Settings
Orbits = OrbNum;
version = 'v1.0';

%Initialize
time = NaN(0,1);
  
  NH = single(NaN(0,1));
  NHcore = single(NaN(0,1));
  NHbeam = single(NaN(0,1));
  NA = single(NaN(0,1));

  VHx = single(NaN(0,1));
  VHy = single(NaN(0,1));
  VHz = single(NaN(0,1));
  VHabs = single(NaN(0,1));

  VHcorex = single(NaN(0,1));
  VHcorey = single(NaN(0,1));
  VHcorez = single(NaN(0,1));

  VAx = single(NaN(0,1));
  VAy = single(NaN(0,1));
  VAz = single(NaN(0,1));
  VAabs = single(NaN(0,1));

  dVHbeamx = single(NaN(0,1));
  dVHbeamy = single(NaN(0,1));
  dVHbeamz = single(NaN(0,1));
  
      TH = single(NaN(0,1)); %Core temperatures
  THperp = single(NaN(0,1));
  THpara = single(NaN(0,1));

      THcore = single(NaN(0,1)); %Total temperatures (treating core and beam as same distribution)
  alphaHcore = single(NaN(0,1));
  THcorePerp = single(NaN(0,1));
  THcorePara = single(NaN(0,1));
  
      TA = single(NaN(0,1));
  alphaA = single(NaN(0,1));
  TAperp = single(NaN(0,1));
  TApara = single(NaN(0,1));

   kappa = single(NaN(0,1));

     Bx = single(NaN(0,1));
     By = single(NaN(0,1));
     Bz = single(NaN(0,1));
   Babs = single(NaN(0,1));

   Valf = single(NaN(0,1)); %Alfven wave speed

   pfit = single(NaN(0,15));

  sigNH = single(NaN(0,1));
  sigNA = single(NaN(0,1));
  sigNHcore = single(NaN(0,1));
  sigNHbeam = single(NaN(0,1)); %NOTE: unit [m^-3] in the V7 files
  
  sigVHx = single(NaN(0,1));
  sigVHy = single(NaN(0,1));
  sigVHz = single(NaN(0,1));

  sigVHcorex = single(NaN(0,1));
  sigVHcorey = single(NaN(0,1));
  sigVHcorez = single(NaN(0,1));

  sigdVHbeamx = single(NaN(0,1));
  sigdVHbeamy = single(NaN(0,1));
  sigdVHbeamz = single(NaN(0,1));

  sigVAx = single(NaN(0,1));
  sigVAy = single(NaN(0,1));
  sigVAz = single(NaN(0,1));

  sigTHcore = single(NaN(0,1));
  sigAlphaHcore = single(NaN(0,1));
  sigAlphaA = single(NaN(0,1));
  sigTA = single(NaN(0,1));

  sigKappa = single(NaN(0,1));

  sigpfit = single(NaN(0,14));
  
  atten_state = single(NaN(0,1));
  isFine = false(0,1);
  fitExitFlag = single(NaN(0,1)); % =2 is good

  frac95 = single(NaN(0,1));

%Main loop
indOrbs = find(ismember(OrbNum,Orbits));
for indNum = 1:length(indOrbs)
    indOrb = indOrbs(indNum);
    

    str = [SWIAswLoc,'mvn_swi_sw_',sprintf('%05d',OrbNum(indOrb)),'_v07_r01.mat'];
    if exist(str,'file')
        SWIAorb = load(str,'time','NH','NHcore','NHbeam','NA','VHx','VHy','VHz','VHabs','VHcorex','VHcorey','VHcorez','VAx','VAy','VAz','VAabs','dVHbeamx','dVHbeamy','dVHbeamz','Valf','THcore','TH','THperp','THpara','TA','alphaHcore','alphaA','THcorePerp','THcorePara','TAperp','TApara','kappa','Bx','By','Bz','Babs','pfit','sigNH','sigNA','sigNHcore','sigNHbeam','sigVHcorex','sigVHcorey','sigVHcorez','sigdVHbeamx','sigdVHbeamy','sigdVHbeamz','sigVHx','sigVHy','sigVHz','sigVAx','sigVAy','sigVAz','sigTHcore','sigTA','sigAlphaHcore','sigAlphaA','sigKappa','isSW','isFine','atten_state','fitExitFlag','frac95');
        
        %Inidices of moments to include
        ind = (SWIAorb.isSW == 1).*(SWIAorb.VHabs > 200) == 1;%.*(SWIAorb.sigNH./SWIAorb.NH < 1) == 1;
        
        time = cat(1,time,SWIAorb.time(ind));
        NH = cat(1,NH,SWIAorb.NH(ind));
        NHcore = cat(1,NHcore,SWIAorb.NHcore(ind));
        NHbeam = cat(1,NHbeam,SWIAorb.NHbeam(ind));
        NA = cat(1,NA,SWIAorb.NA(ind));

        VHx = cat(1,VHx,SWIAorb.VHx(ind));
        VHy = cat(1,VHy,SWIAorb.VHy(ind));
        VHz = cat(1,VHz,SWIAorb.VHz(ind));
        VHabs = cat(1,VHabs,SWIAorb.VHabs(ind));

        VHcorex = cat(1,VHcorex,SWIAorb.VHcorex(ind));
        VHcorey = cat(1,VHcorey,SWIAorb.VHcorey(ind));
        VHcorez = cat(1,VHcorez,SWIAorb.VHcorez(ind));

        VAx = cat(1,VAx,SWIAorb.VAx(ind));
        VAy = cat(1,VAy,SWIAorb.VAy(ind));
        VAz = cat(1,VAz,SWIAorb.VAz(ind));
        VAabs = cat(1,VAabs,SWIAorb.VAabs(ind));

        dVHbeamx = cat(1,dVHbeamx,SWIAorb.dVHbeamx(ind));
        dVHbeamy = cat(1,dVHbeamy,SWIAorb.dVHbeamy(ind));
        dVHbeamz = cat(1,dVHbeamz,SWIAorb.dVHbeamz(ind));

        Valf = cat(1,Valf,SWIAorb.Valf(ind));
        
        TH = cat(1,TH,SWIAorb.TH(ind));
        THperp = cat(1,THperp,SWIAorb.THperp(ind));
        THpara = cat(1,THpara,SWIAorb.THpara(ind));

        THcore = cat(1,THcore,SWIAorb.THcore(ind));
        alphaHcore = cat(1,alphaHcore,SWIAorb.alphaHcore(ind));
        THcorePerp = cat(1,THcorePerp,SWIAorb.THcorePerp(ind));
        THcorePara = cat(1,THcorePara,SWIAorb.THcorePara(ind));
        
        TA = cat(1,TA,SWIAorb.TA(ind));
        alphaA = cat(1,alphaA,SWIAorb.alphaA(ind));
        TAperp = cat(1,TAperp,SWIAorb.TAperp(ind));
        TApara = cat(1,TApara,SWIAorb.TApara(ind));

        kappa = cat(1,kappa,SWIAorb.kappa(ind));
        
        Bx = cat(1,Bx,SWIAorb.Bx(ind));
        By = cat(1,By,SWIAorb.By(ind));
        Bz = cat(1,Bz,SWIAorb.Bz(ind));
        Babs = cat(1,Babs,SWIAorb.Babs(ind));
        
        pfit = cat(1,pfit,SWIAorb.pfit(ind,:));

        sigNH = cat(1,sigNH,SWIAorb.sigNH(ind));
        sigNHcore = cat(1,sigNHcore,SWIAorb.sigNHcore(ind));
        sigNA = cat(1,sigNA,SWIAorb.sigNA(ind));

        sigNHbeam = cat(1,sigNHbeam,SWIAorb.sigNHbeam(ind)); 
        
        sigVHx = cat(1,sigVHx,SWIAorb.sigVHx(ind));
        sigVHy = cat(1,sigVHy,SWIAorb.sigVHy(ind));
        sigVHz = cat(1,sigVHz,SWIAorb.sigVHz(ind));

        sigVHcorex = cat(1,sigVHcorex,SWIAorb.sigVHcorex(ind));
        sigVHcorey = cat(1,sigVHcorey,SWIAorb.sigVHcorey(ind));
        sigVHcorez = cat(1,sigVHcorez,SWIAorb.sigVHcorez(ind));
        
        sigVAx = cat(1,sigVAx,SWIAorb.sigVAx(ind));
        sigVAy = cat(1,sigVAy,SWIAorb.sigVAy(ind));
        sigVAz = cat(1,sigVAz,SWIAorb.sigVAz(ind));

        sigdVHbeamx = cat(1,sigdVHbeamx,SWIAorb.sigdVHbeamx(ind));
        sigdVHbeamy = cat(1,sigdVHbeamy,SWIAorb.sigdVHbeamy(ind));
        sigdVHbeamz = cat(1,sigdVHbeamz,SWIAorb.sigdVHbeamz(ind));

        sigTHcore = cat(1,sigTHcore,SWIAorb.sigTHcore(ind));
        sigAlphaHcore = cat(1,sigAlphaHcore,SWIAorb.sigAlphaHcore(ind));
        sigAlphaA = cat(1,sigAlphaA,SWIAorb.sigAlphaA(ind));

        SWIAorb.sigTA(imag(SWIAorb.sigTA) > 0) = inf;
        SWIAorb.sigTA = real(SWIAorb.sigTA);
        sigTA = cat(1,sigTA,SWIAorb.sigTA(ind));
        
        sigKappa = cat(1,sigKappa,SWIAorb.sigKappa(ind));

        %sigpfit = cat(1,sigpfit,SWIAorb.sigpfit(ind,:)); %Not included in MAVENupstreamV7
        
        atten_state = cat(1,atten_state,SWIAorb.atten_state(ind));
        isFine = cat(1,isFine,SWIAorb.isFine(ind));
        fitExitFlag = cat(1,fitExitFlag,SWIAorb.fitExitFlag(ind));

        frac95 = cat(1,frac95,SWIAorb.frac95(ind));

        disp(['  %% Concatenated: ',str]);
    end
end

VHcoreAbs = sqrt(VHcorex.^2 + VHcorey.^2 + VHcorez.^2);

%Unit conversions
sigNHcore = sigNHcore*1E-6; %[m^-3] -> [cm^-3] %Unit conversion factor was not included in MAVENupstreamV7
sigNHbeam = sigNHbeam*1E-6; %[m^-3] -> [cm^-3] %same

sigVAx = sigVAx/1000;
sigVAy = sigVAy/1000;
sigVAz = sigVAz/1000;

%Fix problem with imaginary uncertainties
sigNH(imag(sigNH) > 0) = sqrt(sigNHcore(imag(sigNH) > 0).^2 + sigNHbeam(imag(sigNH) > 0).^2);

%Normalized magnetic field vector components
Bxn = Bx./Babs;
Byn = By./Babs;
Bzn = Bz./Babs;

%Differential beam velocities along B
dVHbeamB = dVHbeamx.*Bxn + dVHbeamy.*Byn + dVHbeamz.*Bzn;
sigdVHbeamB = sqrt(Bxn.^2.*sigdVHbeamx.^2 + Byn.^2.*sigdVHbeamy.^2 + Bzn.^2.*sigdVHbeamz.^2);

dVHbeamAbs = sqrt(dVHbeamx.^2 + dVHbeamy.^2 + dVHbeamz.^2);

%Calculate differential He2+ velocities
dVAx = VAx-VHcorex;
dVAy = VAy-VHcorey;
dVAz = VAz-VHcorez;

dVAB = dVAx.*Bxn + dVAy.*Byn + dVAz.*Bzn;

sigdVAx = sqrt(sigVAx.^2 + sigVHcorex.^2);
sigdVAy = sqrt(sigVAy.^2 + sigVHcorey.^2);
sigdVAz = sqrt(sigVAz.^2 + sigVHcorez.^2);

sigdVAB = sqrt(Bxn.^2.*sigdVAx.^2 + Byn.^2.*sigdVAy.^2 + Bzn.^2.*sigdVAz.^2);

dVAabs = sqrt(dVAx.^2 + dVAy.^2 + dVAz.^2);
%sigdVAB = sqrt(Bxn.^2.*sigdVAx.^2 + Byn.^2.*sigdVHbeamy.^2 + Bzn.^2.*sigdVHbeamz.^2);

sigTHcorePerp = sqrt((3.*alphaHcore./(1+2*alphaHcore)).^2.^sigTHcore.^2 + (3.*THcore./(1+2*alphaHcore)-(6.*THcore.*alphaHcore)./(1+2*alphaHcore).^2).^2.*sigAlphaHcore.^2);
sigTHcorePara = sqrt((3./(1+2*alphaHcore)).^2.*sigTHcore.^2           + ((-6.*THcore)./(1+2*alphaHcore).^2).^2.*sigAlphaHcore.^2);

sigTAperp = sqrt((3.*alphaA./(1+2*alphaA)).^2.^sigTA.^2 + (3.*TA./(1+2*alphaA)-(6.*TA.*alphaA)./(1+2*alphaA).^2).^2.*sigAlphaA.^2);
sigTApara = sqrt((3./(1+2*alphaA)).^2.*sigTA.^2           + ((-6.*TA)./(1+2*alphaA).^2).^2.*sigAlphaA.^2);

%Elevation of the observed distribution in SWIA's FOV
[azSWIA,elSWIA,~] = cart2sph(pfit(:,4),pfit(:,5),pfit(:,6));
azSWIA = rad2deg(azSWIA); elSWIA = rad2deg(elSWIA);

%Define flags
F1 = ~isFine;
F2 = atten_state == 2;
F3 = frac95 < 0.9;
F4 = ((fitExitFlag < 1.5) + (abs(dVHbeamB)./(Valf/1000) < 0.96) + (kappa > 160) + isnan(sigNHcore) + isnan(sigNHbeam) + isnan(sigNA) + (alphaA < 0.101) + (alphaA > 9.99)) > 0;
F5 = abs(elSWIA) > 30;

%% Save full dataset

timeStr = [datestr(time,'yyyy-mm-ddTHH:MM:SS'),repmat('Z',[length(time) 1])];
dVHbeam = dVHbeamB;
sigdVHbeam = sigdVHbeamB;
dVA = dVAB;
sigdVA = sigdVAB;

sigVHabs = ((sigVHx.^2.*VHx.^2)./(VHx.^2 + VHy.^2 + VHz.^2) + (sigVHy.^2.*VHy.^2)./(VHx.^2 + VHy.^2 + VHz.^2) + (sigVHz.^2.*VHz.^2)./(VHx.^2 + VHy.^2 + VHz.^2)).^(1/2);

%% Save as .mat files
vStr = '01'; %Version number
save([saveLoc,'MVN_SWFITEM_comp_v',vStr,'.mat'],'time','timeStr','NHcore','sigNHcore','NHbeam','sigNHbeam','NA','sigNA','VHcorex','VHcorey','VHcorez','sigVHcorex','sigVHcorey','sigVHcorez','dVHbeam','sigdVHbeam','dVA','sigdVA','THcorePerp','sigTHcorePerp','THcorePara','sigTHcorePara','TAperp','sigTAperp','TApara','sigTApara','kappa','sigKappa','Bx','By','Bz','F1','F2','F3','F4','F5');

save([saveLoc,'MVN_SWFITEM_bulk_v',vStr,'.mat'],'time','timeStr','NH','sigNH','NA','sigNA','VHx','VHy','VHz','sigVHx','sigVHy','sigVHz','VHabs','sigVHabs','TH','TA','sigTA','kappa','sigKappa','Bx','By','Bz','F1','F2','F3','F4','F5');

%% Save as text-files
disp([' %% Saving dataset in ',saveLoc,' as csv text file MVN_SWfit_comp_v',vStr,'.csv']);
disp([' %% Saving dataset in ',saveLoc,' as csv text file MVN_SWfit_bulk_v',vStr,'.csv']);

fid1 = fopen([saveLoc,'MVN_SWFITEM_comp_v',vStr,'.csv'],'w','n','UTF-8');
fid2 = fopen([saveLoc,'MVN_SWFITEM_bulk_v',vStr,'.csv'],'w','n','UTF-8');

fprintf(fid1,'%s\n',['% MVN_SWFITEM_comp_v',vStr,'.csv']); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% This is a dataset of solar wind parameters based on model distribution functions that are fit to phase-space distributions measured in the    '); 
fprintf(fid1,'%s\n','% undisturbed solar wind upstream of Mars by the Solar Wind Ion Analyzer (SWIA), an ion spectrometer onboard the NASA Mars Atmosphere and       '); 
fprintf(fid1,'%s\n','% Volatile EvolutioN (MAVEN) orbiter (Halekas et al. 2015, doi:10.1007/s11214-013-0029-z). Vector quantities are defined in Mars-Sun-Orbit (MSO)'); 
fprintf(fid1,'%s\n','% coordinates where X points to the Sun, Z is perpendicular to the Martian heliocentric orbit plane, Y completes the right-handed system, and   '); 
fprintf(fid1,'%s\n','% the center of the planet defines the origin. For a detailed description of the data reduction method, meanings of parameters, and impacts of  '); 
fprintf(fid1,'%s\n','% flags, see the data descriptor (Ramstad et al. 2025, Scientific Data, doi:TBD).'); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% For use in publications, cite:'); 
fprintf(fid1,'%s\n','%  - The data descriptor: Ramstad et al. (2025), Scientific Data, doi:TBD'); 
fprintf(fid1,'%s\n','%  - The SWIA instrument paper: Halekas et al. (2015), doi:10.1007/s11214-013-0029-z'); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% Variable columns: (name, unit, description)'); 
fprintf(fid1,'%s\n','%  1: UTC-0 time in ISO8601 format (yyyy-mm-ddTHH:MM:SSZ)'); 
fprintf(fid1,'%s\n','%  2:    NHcore [cm^-3], proton (H+) core number density'); 
fprintf(fid1,'%s\n','%  3: sigNHcore [cm^-3], uncertainty on NHcore'); 
fprintf(fid1,'%s\n','%  4:    NHbeam [cm^-3], proton beam number density'); 
fprintf(fid1,'%s\n','%  5: sigNHbeam [cm^-3], uncertainty on NHbeam'); 
fprintf(fid1,'%s\n','%  6:    NA [cm^-3], alpha particle number density (He2+, fully ionized helium)'); 
fprintf(fid1,'%s\n','%  7: sigNA [cm^-3], uncertainty on NA'); 
fprintf(fid1,'%s\n','%  8:    VHcorex [km/s], proton core velocity vector, Cartesian MSO x-component'); 
fprintf(fid1,'%s\n','%  9:    VHcorey [km/s], proton core velocity vector, Cartesian MSO y-component'); 
fprintf(fid1,'%s\n','% 10:    VHcorez [km/s], proton core velocity vector, Cartesian MSO z-component'); 
fprintf(fid1,'%s\n','% 11: sigVHcorex [km/s], uncertainty on VHcorex'); 
fprintf(fid1,'%s\n','% 12: sigVHcorey [km/s], uncertainty on VHcorey'); 
fprintf(fid1,'%s\n','% 13: sigVHcorez [km/s], uncertainty on VHcorez'); 
fprintf(fid1,'%s\n','% 14:    dVHbeam [km/s], differential velocity of the H+ beam along the magnetic field (B) vector, relative to the H+ core velocity'); 
fprintf(fid1,'%s\n','% 15: sigdVHbeam [km/s], uncertainty on dVHbeam'); 
fprintf(fid1,'%s\n','% 16:    dVA [km/s], differential velocity of the alphas along the magnetic field vector, relative to the H+ core velocity'); 
fprintf(fid1,'%s\n','% 17: sigdVA [km/s], uncertainty on dVA'); 
fprintf(fid1,'%s\n','% 18:    THcorePerp [K], perpendicular (relative B) temperature of the core H+ population'); 
fprintf(fid1,'%s\n','% 19: sigTHcorePerp [K], uncertainty on THcorePerp'); 
fprintf(fid1,'%s\n','% 20:    THcorePara [K], parallel (relative B) temperature of the core H+ population'); 
fprintf(fid1,'%s\n','% 21: sigTHcorePara [K], uncertainty on THcorePara'); 
fprintf(fid1,'%s\n','% 22:    TAperp [K], perpendicular (relative B) temperature of the alpha population'); 
fprintf(fid1,'%s\n','% 23: sigTAperp [K], uncertainty on TAperp'); 
fprintf(fid1,'%s\n','% 24:    TApara [K], parallel (relative B) temperature of the alpha population'); 
fprintf(fid1,'%s\n','% 25: sigTApara [K], uncertainty on TAperp'); 
fprintf(fid1,'%s\n','% 26:    kappa [unitless], thermalization factor'); 
fprintf(fid1,'%s\n','% 27: sigKappa [unitless], uncertainty on kappa'); 
fprintf(fid1,'%s\n','% 28:    Bx [nT], magnetic field vector, Cartesian MSO x-component'); 
fprintf(fid1,'%s\n','% 29:    By [nT], magnetic field vector, Cartesian MSO y-component'); 
fprintf(fid1,'%s\n','% 30:    Bz [nT], magnetic field vector, Cartesian MSO z-component'); 
fprintf(fid1,'%s\n','% 31:    F1, Flag 1, =1 if parameters are based on a SWIA coarse distribution'); 
fprintf(fid1,'%s\n','% 32:    F2, Flag 2, =1 if SWIA''s attenuator was in the closed state'); 
fprintf(fid1,'%s\n','% 33:    F3, Flag 3, =1 if likely that the measured distribution was not fit accurately'); 
fprintf(fid1,'%s\n','% 34:    F4, Flag 4, =1 if the fit solution was found at the edge of the constrained parameter space'); 
fprintf(fid1,'%s\n','% 35:    F5, Flag 5, =1 if the solar wind distirbution was observed at an extreme SWIA elevation (abs(el) > 30 deg)'); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% We advice the user to carefully consider the impact of each flag, specified in the data descriptor. '); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% Created by Robin Ramstad (SWIA Deputy PI), Laboratory for Atmospheric and Space Physics, University of Colorado Boulder'); 
fprintf(fid1,'%s\n','%     - for inquiries, contact robin.ramstad@lasp.colorado.edu'); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n',['% This file was created on ',datestr(now,'yyyy-mm-dd')]); 
fprintf(fid1,'%s\n','%'); 
fprintf(fid1,'%s\n','% Variable number 1:,      2:,        3:,      4:,        5:,       6:,      7:,      8:,      9:,     10:,        11:,        12:,        13:,     14:,        15:,     16:,     17:,        18:,           19:,        20:,           21:,     22:,       23:,     24:,       25:,    26:,      27:,      28:,      29:,      30:, 31:, 32:, 33:, 34:, 35:'); 
       strVars1 = '%       Time (UTC-0),  NHcore, sigNHcore,  NHbeam, sigNHbeam,       NA,   sigNA, VHcorex, VHcorey, VHcorez, sigVHcorex, sigVHcorey, sigVHcorez, dVHbeam, sigdVHbeam,     dVA,  sigdVA, THcorePerp, sigTHcorePerp, THcorePara, sigTHcorePara,  TAperp, sigTAperp,  TApara, sigTApara,  kappa, sigKappa,       Bx,       By,       Bz,  F1,  F2,  F3,  F4,  F5';
fprintf(fid1,'%s\n',strVars1); 

fprintf(fid2,'%s\n',['% MVN_FITEM_bulk_v',vStr,'.csv']); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% This is a dataset of bulk solar wind parameters based on model distribution functions that are fit to phase-space distributions measured in   '); 
fprintf(fid2,'%s\n','% the undisturbed solar wind upstream of Mars by the Solar Wind Ion Analyzer (SWIA), an ion spectrometer onboard the NASA Mars Atmosphere and   '); 
fprintf(fid2,'%s\n','% Volatile EvolutioN (MAVEN) orbiter (Halekas et al. 2015, doi:10.1007/s11214-013-0029-z). Vector quantities are defined in Mars-Sun-Orbit (MSO)'); 
fprintf(fid2,'%s\n','% coordinates where X points to the Sun, Z is perpendicular to the Martian heliocentric orbit plane, Y completes the right-handed system, and   '); 
fprintf(fid2,'%s\n','% the center of the planet defines the origin. For a detailed description of the data reduction method, meanings of parameters, and impacts of  '); 
fprintf(fid2,'%s\n','% flags, see the data descriptor (Ramstad et al. 2025, Scientific Data, doi:TBD).'); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% For use in publications, cite:'); 
fprintf(fid2,'%s\n','%  - The data descriptor: Ramstad et al. (2025), Scientific Data, doi:TBD'); 
fprintf(fid2,'%s\n','%  - The SWIA instrument paper: Halekas et al. (2015), doi:10.1007/s11214-013-0029-z'); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% Variable columns: (name, unit, description)'); 
fprintf(fid2,'%s\n','%  1: UTC-0 time in ISO8601 format (yyyy-mm-ddTHH:MM:SSZ)'); 
fprintf(fid2,'%s\n','%  2:    NH [cm^-3], proton (H+) number density'); 
fprintf(fid2,'%s\n','%  3: sigNH [cm^-3], uncertainty on NH'); 
fprintf(fid2,'%s\n','%  4:    NA [cm^-3], alpha particle number density (He2+, fully ionized helium)'); 
fprintf(fid2,'%s\n','%  5: sigNA [cm^-3], uncertainty on NA'); 
fprintf(fid2,'%s\n','%  6:    VHx [km/s], proton bulk velocity vector, Cartesian MSO x-component'); 
fprintf(fid2,'%s\n','%  7:    VHy [km/s], proton bulk velocity vector, Cartesian MSO y-component'); 
fprintf(fid2,'%s\n','%  8:    VHz [km/s], proton bulk velocity vector, Cartesian MSO z-component'); 
fprintf(fid2,'%s\n','%  9: sigVHx [km/s], uncertainty on VHx'); 
fprintf(fid2,'%s\n','% 10: sigVHy [km/s], uncertainty on VHy'); 
fprintf(fid2,'%s\n','% 11: sigVHz [km/s], uncertainty on VHz'); 
fprintf(fid2,'%s\n','% 12:  VHabs [km/s], bulk proton speeed'); 
fprintf(fid2,'%s\n','% 13: sigVHabs [km/s], uncertainty on dVHbeam'); 
fprintf(fid2,'%s\n','% 14:    TH [K], scalar temperature of the bulk H+ distribution'); 
fprintf(fid2,'%s\n','% 15:    TA [K], scalar temperature of the bulk He2+ distribution'); 
fprintf(fid2,'%s\n','% 16: sigTA [K], uncertainty on TA'); 
fprintf(fid2,'%s\n','% 17:    kappa [unitless], thermalization factor'); 
fprintf(fid2,'%s\n','% 18: sigKappa [unitless], uncertainty on kappa'); 
fprintf(fid2,'%s\n','% 19:    Bx [nT], magnetic field vector, Cartesian MSO x-component'); 
fprintf(fid2,'%s\n','% 20:    By [nT], magnetic field vector, Cartesian MSO y-component'); 
fprintf(fid2,'%s\n','% 21:    Bz [nT], magnetic field vector, Cartesian MSO z-component'); 
fprintf(fid2,'%s\n','% 22:    F1, Flag 1, =1 if parameters are based on a SWIA coarse distribution'); 
fprintf(fid2,'%s\n','% 23:    F2, Flag 2, =1 if SWIA''s attenuator was in the closed state'); 
fprintf(fid2,'%s\n','% 24:    F3, Flag 3, =1 if likely that the measured distribution was not fit accurately'); 
fprintf(fid2,'%s\n','% 25:    F4, Flag 4, =1 if the fit solution was found at the edge of the constrained parameter space'); 
fprintf(fid2,'%s\n','% 26:    F5, Flag 5, =1 if the solar wind distirbution was observed at an extreme SWIA elevation (abs(el) > 30 deg)'); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% We advice the user to carefully consider the impact of each flag, specified in the data descriptor. '); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% Created by Robin Ramstad (SWIA Deputy PI), Laboratory for Atmospheric and Space Physics, University of Colorado Boulder'); 
fprintf(fid2,'%s\n','%     - for inquiries, contact robin.ramstad@lasp.colorado.edu'); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n',['% This file was created on ',datestr(now,'yyyy-mm-dd')]); 
fprintf(fid2,'%s\n','%'); 
fprintf(fid2,'%s\n','% Variable number 1:,      2:,      3:,      4:,      5:,      6:,      7:,      8:,      9:,     10:,     11:,     12:,      13:,     14:,     15:,     16:,     17:,      18:,      19:,      20:,      21:, 22:, 23:, 24:, 25:, 26:'); 
       strVars2 = '%       Time (UTC-0),      NH,   sigNH,      NA,   sigNA,     VHx,     VHy,     VHz,  sigVHx,  sigVHy,  sigVHz,   VHabs, sigVHabs,      TH,      TA,   sigTA,   kappa, sigKappa,       Bx,       By,       Bz,  F1,  F2,  F3,  F4,  F5';
fprintf(fid2,'%s\n',strVars2); 

dispThres = 0; tic; %Progress update threshold
for ii = 1:length(time)
    
    if ii/length(time) >= dispThres
        dispThres = dispThres+0.05;
        disp(['   %% Exporting to csv text files: ',sprintf('%.0f',100*ii/length(time)),'% completed (',sprintf('%.1f',toc/60),' min)']);
    end

    %Strings for exponential notation of variables, remove unnecessary zeros in exponent to save space
    strNHcore     = sprintf('%7.2E',   NHcore(ii)); if ~isnan(NHcore(ii)); strNHcore(end-1) = []; end
    strSigNHcore  = sprintf('%7.2E',sigNHcore(ii)); if ~isnan(sigNHcore(ii)) && ~isinf(sigNHcore(ii)); strSigNHcore(end-1) = []; end
    strNHbeam     = sprintf('%7.2E',   NHbeam(ii)); if ~isnan(NHbeam(ii)); strNHbeam(end-1) = []; end
    strSigNHbeam  = sprintf('%7.2E',sigNHbeam(ii)); if ~isnan(sigNHbeam(ii)) && ~isinf(sigNHbeam(ii)); strSigNHbeam(end-1) = []; end
    strNA         = sprintf('%7.2E',   NA(ii)); if ~isnan(NA(ii)) && ~isinf(NA(ii)); strNA(end-1) = []; end
    strSigNA      = sprintf('%7.2E',sigNA(ii)); if ~isnan(sigNA(ii)) && ~isinf(sigNA(ii)); strSigNA(end-1) = []; end
    strVHcorex    = sprintf('%+7.1f',   VHcorex(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHcorex = sprintf('%7.1f',sigVHcorex(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strVHcorey    = sprintf('%+7.1f',   VHcorey(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHcorey = sprintf('%7.1f',sigVHcorey(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strVHcorez    = sprintf('%+7.1f',   VHcorez(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHcorez = sprintf('%7.1f',sigVHcorez(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strdVHbeam    = sprintf('%+7.1f',   dVHbeam(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigdVHbeam = sprintf('%7.1f', sigdVHbeam(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strdVA        = sprintf('%+7.1f',   dVA(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigdVA     = sprintf('%7.1f', sigdVA(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strTH    = sprintf('%7.2E',TH(ii)); if ~isnan(TH(ii)); strTH(end-1) = []; end
    strTHcorePerp    = sprintf('%7.2E',THcorePerp(ii)); if ~isnan(THcorePerp(ii)); strTHcorePerp(end-1) = []; end
    strSigTHcorePerp = sprintf('%7.2E',sigTHcorePerp(ii)); if ~isnan(sigTHcorePerp(ii)) && ~isinf(sigTHcorePerp(ii)); strSigTHcorePerp(end-1) = []; end
    strTHcorePara    = sprintf('%7.2E',THcorePara(ii)); if ~isnan(THcorePara(ii)); strTHcorePara(end-1) = []; end
    strSigTHcorePara = sprintf('%7.2E',sigTHcorePara(ii)); if ~isnan(sigTHcorePara(ii)) && ~isinf(sigTHcorePara(ii)); strSigTHcorePara(end-1) = []; end
    strTAperp     = sprintf('%7.2E',TAperp(ii)); if ~isnan(TAperp(ii)); strTAperp(end-1) = []; end
    strSigTAperp  = sprintf('%7.2E',sigTAperp(ii)); if ~isnan(sigTAperp(ii)) && ~isinf(sigTAperp(ii)); strSigTAperp(end-1) = []; end
    strTApara     = sprintf('%7.2E',TApara(ii)); if ~isnan(TApara(ii)); strTApara(end-1) = []; end
    strSigTApara  = sprintf('%7.2E',sigTApara(ii)); if ~isnan(sigTApara(ii)) && ~isinf(sigTApara(ii)); strSigTApara(end-1) = []; end
    strKappa      = sprintf('%7.2f',kappa(ii));%if ~isnan(kappa(ii)); kappa(end-1) = []; end
    strSigKappa   = sprintf('%7.2E',sigKappa(ii)); if ~isnan(sigKappa(ii)) && ~isinf(sigKappa(ii)); strSigKappa(end-1) = []; end
    strBx         = sprintf('%+8.2E',Bx(ii)); if ~isnan(Bx(ii)); strBx(end-1) = []; end
    strBy         = sprintf('%+8.2E',By(ii)); if ~isnan(By(ii)); strBy(end-1) = []; end
    strBz         = sprintf('%+8.2E',Bz(ii)); if ~isnan(Bz(ii)); strBz(end-1) = []; end
    
    %Cap any extremely large uncertainties
    if length(strSigNHcore) > 7;  strSigNHcore  = '    Inf'; end
    if length(strSigNHbeam) > 7;  strSigNHbeam  = '    Inf'; end
    if length(strSigNA) > 7;      strSigNA      = '    Inf'; end
    if length(strSigVHcorex) > 7; strSigVHcorex = '    Inf'; end
    if length(strSigVHcorey) > 7; strSigVHcorey = '    Inf'; end
    if length(strSigVHcorez) > 7; strSigVHcorez = '    Inf'; end
    if length(strSigdVHbeam) > 7; strSigdVHbeam = '    Inf'; end
    if length(strSigdVA) > 7;     strSigdVA     = '    Inf'; end
    if length(strSigTHcorePerp) > 7;  strSigTHcorePerp = '    Inf'; end
    if length(strSigTHcorePara) > 7;  strSigTHcorePara = '    Inf'; end
    if length(strSigTAperp) > 7;  strSigTAperp = '    Inf'; end
    if length(strSigTApara) > 7;  strSigTApara = '    Inf'; end
    if length(strSigKappa) > 7;   strSigKappa  = '    Inf'; end

    strNH       = sprintf('%7.2E',   NH(ii)); if ~isnan(NH(ii)); strNH(end-1) = []; end
    strSigNH    = sprintf('%7.2E',sigNH(ii)); if ~isnan(sigNH(ii)) && ~isinf(sigNH(ii)); strSigNH(end-1) = []; end
    strVHx    = sprintf('%+7.1f',   VHx(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHx = sprintf('%7.1f',sigVHx(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strVHy    = sprintf('%+7.1f',   VHy(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHy = sprintf('%7.1f',sigVHy(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strVHz    = sprintf('%+7.1f',   VHz(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strSigVHz = sprintf('%7.1f',sigVHz(ii)); %if ~isnan(VHcorex(ii)); strVHcorex(end-1) = []; end
    strVHabs    = sprintf('%7.1f',VHabs(ii));
    strSigVHabs = sprintf('%8.1f',sigVHabs(ii));
    strTH       = sprintf('%7.2E',   TH(ii)); if ~isnan(TH(ii)) && ~isinf(TH(ii)); strTH(end-1) = []; end
    strTA       = sprintf('%7.2E',   TA(ii)); if ~isnan(TA(ii)) && ~isinf(TA(ii)); strTA(end-1) = []; end
    strSigTA    = sprintf('%7.2E',sigTA(ii)); if ~isnan(sigTA(ii)) && ~isinf(sigTA(ii)); strSigTA(end-1) = []; end
    
    if length(strSigNH) > 7;    strSigNH  = '    Inf'; end
    if length(strSigVHx) > 7;   strSigVHx = '    Inf'; end
    if length(strSigVHy) > 7;   strSigVHy = '    Inf'; end
    if length(strSigVHz) > 7;   strSigVHz = '    Inf'; end
    if length(strSigVHabs) > 8; strSigVHabs = '     Inf'; end

    %Print to full components file
    printStr1 = [timeStr(ii,:),', ',strNHcore,',   ',strSigNHcore,', ',strNHbeam,',    ',strSigNHbeam,', ',strNA,', ',strSigNA,', ',strVHcorex,', ',strVHcorey,', ',strVHcorez,',    ',strSigVHcorex,',    ',strSigVHcorey,',    ',strSigVHcorez,', ',strdVHbeam,',    ',strSigdVHbeam,', ',strdVA,', ',strSigdVA,',    ',strTHcorePerp,',       ',strSigTHcorePerp,',    ',strTHcorePara,',       ',strSigTHcorePara,', ',strTAperp,',   ',strSigTAperp,', ',strTApara,',   ',strSigTApara,',',strKappa,',  ',strSigKappa,', ',strBx,', ',strBy,', ',strBz,', ',sprintf('%3.0f',F1(ii)),', ',sprintf('%3.0f',F2(ii)),', ',sprintf('%3.0f',F3(ii)),', ',sprintf('%3.0f',F4(ii)),', ',sprintf('%3.0f',F5(ii))];
    
    fprintf(fid1,'%s\n',printStr1); %Print each line (time-stamp of moments)
    
    %Print to bulk parameter file
    printStr2 = [timeStr(ii,:),', ',strNH,', ',strSigNH,', ',strNA,', ',strSigNA,', ',strVHx,', ',strVHy,', ',strVHz,', ',strSigVHx,', ',strSigVHy,', ',strSigVHz,', ',strVHabs,', ',strSigVHabs,', ',strTH,', ',strTA,', ',strSigTA,', ',strKappa,',  ',strSigKappa,', ',strBx,', ',strBy,', ',strBz,', ',sprintf('%3.0f',F1(ii)),', ',sprintf('%3.0f',F2(ii)),', ',sprintf('%3.0f',F3(ii)),', ',sprintf('%3.0f',F4(ii)),', ',sprintf('%3.0f',F5(ii))];
    fprintf(fid2,'%s\n',printStr2); %Print each line (time-stamp of moments)

    if (length(printStr1) ~= length(strVars1)) || (length(printStr2) ~= length(strVars2))
        warning(['  %% Export to .csv file FAILED (ii = ',num2str(ii),')']);
        fclose(fid1);
        fclose(fid2);
        return;
    end
    
end
fclose(fid1);
fclose(fid2);
