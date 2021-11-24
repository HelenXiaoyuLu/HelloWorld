%% Collect all power correction table 
% As of 09/10/2020
path = 'C:\Users\13801\OneDrive\Documents\Thorlabs\Optical Power Monitor';
fname = '';
T1 = load('tra mGold\20200215_peakvaluest_jobs_700-20-1080');
T = table();
T.Wavelength = T1.peakvaluest(1,:)';
T.("Power at Sample Plane") = T1.peakvaluest(2,:)';
T_power_20200215_jobs_20nm = T;

%%
T8 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200215 Spectra mGold\20200215_peakvaluest_manual_700-20-1080');
T = table();
T.Wavelength = T8.peakvaluest(1,:)';
T.("Power at Sample Plane") = T8.peakvaluest(2,:)';
T_power_20200215_manual_20nm_10mW = T;

%%
T2 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200227 Spectra BFP LP out\20200227_JOBs_peakvaluest.mat');
T = table();
T.Wavelength = T2.peakvaluest(1,:)';
T.("Power at Sample Plane") = T2.peakvaluest(2,:)';
T_power_20200227_jobs_20nm = T;

%%
T3 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200305 Spectra Final\20200305_peakvaluest_jobs_4_100mW');
T = table();
T.Wavelength = T3.peakvaluest(1,:)';
T.("Power at Sample Plane") = T3.peakvaluest(2,:)';
T_power_20200305_jobs_20nm_100mW = T;

%%
T4 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200305 Spectra Final\20200305_peakvaluest_jobs_700-10-1080');
T = table();
T.Wavelength = T4.peakvaluest(1,:)';
T.("Power at Sample Plane") = T4.peakvaluest(2,:)';
T_power_20200305_jobs_10nm_10mW = T;

%%
T5 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200305 Spectra Final\20200305_peakvaluest_manual_700-10-1080');
T = table();
T.Wavelength = T5.peakvaluest(1,:)';
T.("Power at Sample Plane") = T5.peakvaluest(2,:)';
T_power_20200305_manual_10nm_10mW = T;

%%
T6 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200904 ramping power\20200904_peakvaluest_jobs_10mwrange');
T = table();
T.Wavelength = T6.peakvaluest(1,:)';
T.("Power at Sample Plane") = T6.peakvaluest(2,:)';
T_power_20200904_manual_10nm_10mW = T;

%%
T7 = load('D:\OneDrive - Baylor College of Medicine\Paper_201906_GEVI\Spectra\20200904 ramping power\20200904_peakvaluest_jobs_100mwrange');
T = table();
T.Wavelength = T7.peakvaluest(1,:)';
T.("Power at Sample Plane") = T7.peakvaluest(2,:)';
T_power_20200904_manual_10nm_100mW = T;

%%
save('D:\Github\FieldStim\Visual\+FOVbase_2PSpectra\T_power_collection.mat',...
 'T_power_20200215_jobs_20nm','T_power_20200215_manual_20nm_10mW',...
 'T_power_20200227_jobs_20nm','T_power_20200305_jobs_10nm_10mW',...
 'T_power_20200305_jobs_20nm_100mW','T_power_20200305_manual_10nm_10mW',...
 'T_power_20200904_manual_10nm_100mW','T_power_20200904_manual_10nm_10mW');
