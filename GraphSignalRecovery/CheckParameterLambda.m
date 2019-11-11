
function lambda = CheckParameterLambda(param,BB,L,BC)

    lambdacheck = fullfile(param.EigPath, ['RegularizationParameters_LCurve_',param.subject,'_',num2str(param.smoothR),'.mat']);
    
    if ~exist(lambdacheck, 'file')
        
        % Load preliminary fMRI volume
        Inputs_fMRI
        
        i = 20; % frame to evaluate
        vol = V(i,:);
        
        clear V
        
        disp('Calculating correct parameter lambda..')
        lambdas = 1:60;
%         diffs = zeros(length(lambdas),1);
        x1norms = zeros(length(lambdas),1);
        parfor c = 1:length(lambdas)
            lambda = lambdas(c);
            L_mat = BB+lambda*L;
            L1 = ichol(L_mat);
            xsignal = zeros(size(L1,1),1);
            xsignal(BC) = vol(BC);
            [x1,~,~,~,~] = pcg(L_mat,xsignal,1e-8,100,L1,L1');
            x1norms(c) = norm(x1'*(L*x1));
%             diff = B*x1-xsignal;
%             diffs(c) = norm(diff)^2;
        end

        lambdaStars.lambdas = lambdas;
        lambdaStars.x1norms = x1norms;
        save(fullfile(param.EigPath, ['RegularizationParameters_LCurve_',param.subject,'.mat']),'lambdaStars')
        
    else
        load(lambdacheck)
        x1norms = lambdaStars.x1norms;
        lambdas = lambdaStars.lambdas;
    end
    
    lambda = knee_pt(x1norms,lambdas);
end