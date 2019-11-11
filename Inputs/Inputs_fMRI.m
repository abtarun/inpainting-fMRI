%%%%%% INPUT FUNCTIONAL VOLUMES PATH

% filepath of preprocessed functional volumes:
% should already be in matrix form (timepoints x number of voxel)

param.functional = fullfile(param.HCPDatapath, param.subject,'functional');

if ~isfield(param,'session')
    funcName = 'BCVolumes_s6_100307_session_rfMRI_REST1_LR.mat';
else
    funcName = ['BCVolumes_s6_100307_session_',param.session,'.mat'];
end

% to detrend volumes (put 1 if volumes were not detrended in preprocessing)
Detrend  = 0;

prefix = 's6';
% Contains all inputs regarding functional MRI volumes which are subject-specific
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Functional volumes 

% Loading functional data matrix
load(fullfile(param.functional, funcName))

directory = dir(fullfile(param.functional, [prefix,'*']));
fHeader = spm_vol(fullfile(param.functional,directory(1).name));
hdr=cbiReadNiftiHeader(fHeader.fname);


% If the functional volumes were not detrended, we perform detrending
% signal inpainting
if Detrend
    disp('Let us detrend volume first..')
    CUTNUMBER=10;
    SegmentLength = ceil(size(V,2) / CUTNUMBER);
    for iCut=1:CUTNUMBER
        if iCut~=CUTNUMBER
            Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
        else
            Segment = (iCut-1)*SegmentLength+1 : size(V,2);
        end
        V(:,Segment) = detrend(V(:,Segment));
        V(:,Segment) = detrend(V(:,Segment),'constant');
        fprintf('.');
    end
end
NumScans = size(V,1);