
% Contains the input parameters for the construction of the brain grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specifies the parameters of the braingraph and the eigendecomposition
param.odfPow = 40; % ODF power.
param.neighborhood = 3; % 3 or 5 [This is for white matter; gray matter always 3]
param.N_fibers = 5; % dsi studio's and M&L's default: 5
% Tune FA/QA
param.alpha = 2;

% Choose the parameters for the eigendecomposition of the Laplacian to the 
% obtain eigenmodes..
% whether to have a percent bandwidth or constant
param.percent = 0;
param.bandwidth = 1e-3; 
param.c_bandwidth = 5; % just to visualize first 5 eigenmodes
param.opts.issym=1;
param.opts.isreal=1;
param.opts.maxit=2500;
param.opts.disp=1;
param.normalize = 1;
param.normalize_type = 1;


% Specifies where to save all results of the brain grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.title_ODF = ['ODF_Neigh_',num2str(param.neighborhood),'_ODFPower_',num2str(param.odfPow)];
param.SaveDirectory = fullfile(param.HCPDatapath, 'BrainGraph_results', param.subject,param.title_ODF);
param.structural = fullfile(param.HCPDatapath, param.subject, 'T1w');
if ~exist(param.SaveDirectory, 'dir')
    mkdir(param.SaveDirectory)
end