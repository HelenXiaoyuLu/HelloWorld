%% Overwrite offline access folder to database
   gc = util.globalConfig.getInstance();
   gc.removeConfig("");
   gc.setConfig('User', ["Asuka", "Harry", "James", "Eric", "Xiaoyu", "Kevin", "Amir", "Guest"]);
   gc.setConfig('Database\server', "127.0.0.1");
   gc.setConfig('Database\user', "root");%   gc.setConfig('Database\pwd', "voltage");
   gc.setConfig('Version', "20.0601");
   gc.setConfig('Path\DJ_DB', 'E:\MongoDB\cache');
   gc.setConfig('Path\ARCHIVE', '\\STPIERRELAB7910\E');
   gc.setConfig('Path\TEST', '\\STPIERRELAB7910\E\Images\DEMO');
   gc.setConfig('Path\DB_OFFLINE', 'D:\OneDrive - Baylor College of Medicine\MATLABDB')
   gc.commit();
   
%% Leaving debug mode   
gc.loadConfig('config.xml')
