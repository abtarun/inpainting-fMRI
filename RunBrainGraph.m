function param = RunBrainGraph(param)


    % Date and time when the routines are called
    param.date = strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
    
    param.constW = param.c_bandwidth;
    % Original parameters
    param_init = param;
    
   
        
    disp(['Extracting whole brain graph using ODF-based design for ..',param.subject])


    % Calls the mdh routine to construct the graph and return adjacency matrix

    [param.WB, param.G] = main_braingraph(param);

    [param.WB.A,param.WB.D]=slepNormalize(param.WB.A,param.normalize, param.normalize_type);

    % Computing the eigendecomposition of the Laplacian

    if param.percent
        param.constW = round(param.bandwidth*param.G.N_wb);
    else
        param.constW = param.c_bandwidth;
    end

    disp('Taking the eigenvalues...')

    [param.WB.Utr,param.WB.S]=slepEigsLaplacian(param.WB.A,param.constW,param.opts);

    % Saves the eigenvalues and eigenvectors

    SaveValues(param)

    SaveToNifti(param, param.WB.Utr,'Spectrum_WB')


end