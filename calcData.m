%Herl-Keeogo Overground Protocol
%Author: Maddi Viteri
%Last Updated: 02/21/2025 (MATLAB R2024a)

%this program does the calculations for the ASB 2025 abstract
%medain and IQR for APDisp (mm), MLDisp(mm), Vel(MM),
%time2_maxVel (% for total foot)

%-------------------------------------------------------------------------%

function [disp, Vel, t2_maxVel] = calcData(noraxon,trialNames)

    max_APDisp = 0; min_APDisp = 0;
    max_MLDisp = 0; min_MLDisp = 0;
    median_Vel = 0; iqr_Vel = 0;
    median_t2MaxVel = 0; iqr_t2MaxVel = 0;

    %extract and format data

    NONSS01 = noraxon(1); NONSS02 = noraxon(2);
    NONFS01 = noraxon(3); NONFS02 = noraxon(4);
    KONSS01 = noraxon(5); KONSS02 = noraxon(6);
    KONFS01 = noraxon(7); KONFS02 = noraxon(8);

    NonSS01L = NONSS01{1,1}{1,3};NonSS01R = NONSS01{1,1}{1,4};
    NonSS02L = NONSS02{1,1}{1,3};NonSS02R = NONSS02{1,1}{1,4};
    NonFS01L = NONFS01{1,1}{1,3};NonFS01R = NONFS01{1,1}{1,4};
    NonFS02L = NONFS02{1,1}{1,3};NonFS02R = NONFS02{1,1}{1,4};

    KonSS01L = KONSS01{1,1}{1,3};KonSS01R = KONSS01{1,1}{1,4};
    KonSS02L = KONSS02{1,1}{1,3};KonSS02R = KONSS02{1,1}{1,4};
    KonFS01L = KONFS01{1,1}{1,3};KonFS01R = KONFS01{1,1}{1,4};
    KonFS02L = KONFS02{1,1}{1,3};KonFS02R = NONFS02{1,1}{1,4};

    %displacement calcs

    NonSSLX = vertcat(NonSS01L.COP_X,NonSS02L.COP_X); NonSSLY = vertcat(NonSS01L.COP_Y,NonSS02L.COP_Y);
    NonFSLX = vertcat(NonFS01L.COP_X,NonFS02L.COP_X); NonFSLY = vertcat(NonFS01L.COP_Y,NonFS02L.COP_Y);
    KonSSLX = vertcat(KonSS01L.COP_X,KonSS02L.COP_X); KonSSLY = vertcat(KonSS01L.COP_Y,KonSS02L.COP_Y);
    KonFSLX = vertcat(KonFS01L.COP_X,KonFS02L.COP_X); KonFSLY = vertcat(KonFS01L.COP_Y,KonFS02L.COP_Y);

    f_NonSSLX = unique(NonSSLX); f_NonSSLY = unique(NonSSLY);
    f_NonFSLX = unique(NonFSLX); f_NonFSLY = unique(NonFSLY);
    f_KonSSLX = unique(KonSSLX); f_KonSSLY = unique(KonSSLY);
    f_KonFSLX = unique(KonFSLX); f_KonFSLY = unique(KonFSLY);

    NonSS_Disp = struct('maxAP',max(f_NonSSLY),'minAP',min(f_NonSSLY),'diffAP',max(f_NonSSLY)-min(f_NonSSLY),'maxML',max(f_NonSSLX),'minML',min(f_NonSSLX),'diffML',max(f_NonSSLX)-min(f_NonSSLX));
    NonFS_Disp = struct('maxAP',max(f_NonFSLY),'minAP',min(f_NonFSLY),'diffAP',max(f_NonFSLY)-min(f_NonFSLY),'maxML',max(f_NonFSLX),'minML',min(f_NonFSLX),'diffML',max(f_NonFSLX)-min(f_NonFSLX));
    KonSS_Disp = struct('maxAP',max(f_KonSSLY),'minAP',min(f_KonSSLY),'diffAP',max(f_KonSSLY)-min(f_KonSSLY),'maxML',max(f_KonSSLX),'minML',min(f_KonSSLX),'diffML',max(f_KonSSLX)-min(f_KonSSLX));
    KonFS_Disp = struct('maxAP',max(f_KonFSLY),'minAP',min(f_KonFSLY),'diffAP',max(f_KonFSLY)-min(f_KonFSLY),'maxML',max(f_KonFSLX),'minML',min(f_KonFSLX),'diffML',max(f_KonFSLX)-min(f_KonFSLX));

    disp = [NonSS_Disp,NonFS_Disp,KonSS_Disp,KonFS_Disp]; %mm

    %velocity calculations
    
    Nss01L_vel = (diff(NonSS01L.COP_X)./1000)./diff(NonSS01L.time); %mm to m
    Nss01R_vel = (diff(NonSS01R.COP_Y)./1000)./diff(NonSS01R.time); %mm to m
    Nss01_vel = sqrt(Nss01L_vel.^2 + Nss01R_vel.^2);
    Nss01_Mvel = median(Nss01_vel,"omitnan"); 
    Nss01_iqrV = iqr(Nss01_vel);
    

    Nss02L_vel = (diff(NonSS02L.COP_X)./1000)./diff(NonSS02L.time); %mm to m
    Nss02R_vel = (diff(NonSS02R.COP_Y)./1000)./diff(NonSS02R.time); %mm to m
    Nss02_vel = sqrt(Nss02L_vel.^2 + Nss02R_vel.^2);
    Nss02_Mvel = median(Nss02_vel,"omitnan"); 
    Nss02_iqrV = iqr(Nss02_vel);

    timeArrays = {NonSS01L.time(1:end-1, :),NonSS02L.time(1:end-1, :)};
    velocityArrays = {Nss01_vel,Nss02_vel};

    totalmedian_velocity = calcMedianVelocity(timeArrays, velocityArrays);




    Kss01L_vel = (diff(KonSS01L.COP_X)./1000)./diff(KonSS01L.time); %mm to m
    Kss01R_vel = (diff(KonSS01R.COP_Y)./1000)./diff(KonSS01R.time); %mm to m
    Kss01_vel = sqrt(Kss01L_vel.^2 + Kss01R_vel.^2);
    Kss01_Mvel = median(Kss01_vel,"omitnan"); 
    Kss01_iqrV = iqr(Kss01_vel);

    Kss02L_vel = (diff(KonSS02L.COP_X)./1000)./diff(KonSS02L.time); %mm to m
    Kss02R_vel = (diff(KonSS02R.COP_Y)./1000)./diff(KonSS02R.time); %mm to m
    Kss02_vel = sqrt(Kss02L_vel.^2 + Kss02R_vel.^2);
    Kss02_Mvel = median(Kss02_vel,"omitnan"); 
    Kss02_iqrV = iqr(Kss02_vel);


    NonSS_vel = struct('median',Nss01_Mvel,'IQR',Nss01_iqrV);
    NonFS_vel = struct('median',Nss02_Mvel,'IQR',Nss02_iqrV);
    KonSS_vel = struct('median',Kss01_Mvel,'IQR',Kss01_iqrV);
    KonFS_vel = struct('median',Kss02_Mvel,'IQR',Kss02_iqrV);

    Vel = [NonSS_vel, NonFS_vel, KonSS_vel, KonFS_vel];

    
    %time 2 max velocity

    t2_maxVel = struct('median',median_t2MaxVel,'IQR',iqr_t2MaxVel);

end%end function