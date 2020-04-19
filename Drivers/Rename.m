% Pure Rename 
OriginalName = 'HEK_JEDI-2P_26-*.nd2';
ParentFolder = 'D:\OneDrive - rice.edu\Francois\Paper\Data\Power measurements\PowerSpectrum\20200227 Spectrum\1080-700-1080 BFP LP out\GEVIs';
dirtable = dir(fullfile(ParentFolder,OriginalName));
NewName = 'Kir_JEDI-2P';
for j = 1:length(dirtable)
    [oldPath, oldfName, fext] = fileparts(strcat(dirtable(j).folder,'\',dirtable(j).name));
    varName = split(oldfName,'_');
    newfName = strcat(NewName,'_', varName{3});
    movefile(fullfile(oldPath, [oldfName, fext]), fullfile(oldPath, [newfName, fext]), 'f');
end