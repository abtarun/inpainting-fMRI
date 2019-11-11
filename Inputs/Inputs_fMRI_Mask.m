
% Paths

param.EigPath = fullfile(param.HCPDatapath, param.subject, 'T1w','graph');

param.SaveInpainting = fullfile(param.HCPDatapath, 'Inpainting_results');

if ~exist(param.SaveInpainting,'dir')
    mkdir(param.SaveInpainting)
end

% Reads gray matter mask from segmented T1 image
Mask = spm_read_vols(spm_vol(fullfile(param.structural,...
    'c1T1w_acpc_dc_restore_1.25.nii')));

Mask(Mask<0.3)=0;

Mask = bwareaopen(logical(Mask),30);

Mask = imfill(Mask,'holes');

BCBrainMask = logical(Mask);

param.smoothR = 6;

load(fullfile(param.EigPath, 'indices_wb.mat'))

param.indices = indices_wb;

BC = find(BCBrainMask);

[sim_ind,~,BC] = intersect(BC, param.indices);

if ~isfield(param.WB,'A')
    load(fullfile(param.EigPath,'A_wb.mat')) 
end


