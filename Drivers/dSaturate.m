%% discard saturated pixels, match histogram and save as tiff --> later masking
% Set path: parental and subfolders
% ImgsDir: list of all imgs to be processed and converted
Parent = 'E:\Images\Xiaoyu\20200708 photoactivation train';
subJobs = '\TrainStim 5pctGreen';
subfNames = '\*.nd2';
ImgsDir = dir(strcat(Parent,subJobs,subfNames));
Helpers.stackdSat(ImgsDir,true)
