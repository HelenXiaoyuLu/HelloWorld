%% discard saturated pixels, match histogram and save as tiff --> later masking
% Set path: parental and subfolders
% ImgsDir: list of all imgs to be processed and converted
Parent = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\';
subJobs = 'HV20';
subfNames = '\*.nd2';
ImgsDir = dir(strcat(Parent,subJobs,subfNames));
stackdSat(ImgsDir,true)
