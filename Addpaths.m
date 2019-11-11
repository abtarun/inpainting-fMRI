function param = Addpaths
    % Path to the spm12
    addpath(genpath('/Users/anjalitarun/Documents/SPM/spm12'));
    
    addpath(genpath(pwd));
    
    
    % Code base path
    param.DSIstudio = '/Users/anjalitarun/Documents/Softwares/DSI_studio/Contents/MacOS';

    % Database path
    param.HCPDatapath = '/Users/anjalitarun/Downloads/HCP';
    
    
end