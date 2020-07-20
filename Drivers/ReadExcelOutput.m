%%
clear
clc
dirpath ='D:\OneDrive - rice.edu\Francois\RedGEVIs\Wet\Results\XLRED001-005\201907*1000Hz.xlsx';
controlName = ['ASAP1','JEDI-1P','JEDI-1P-VSD-N-1-C-1-cpmOrange2-174-173'];
libStore = struct;
libFiles = dir(dirpath);
for i = 1:length(libFiles)
    libStore(i).table = readtable(fullfile(libFiles(i).folder,libFiles(i).name),'Sheet','Summary (Library 1)');
end

