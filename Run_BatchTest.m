% Run a batch of tests and analyse performance.

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% How many tests?
num_tests = 20;

Stats = struct('MMSE', zeros(num_tests,1), 'LostTracks', zeros(num_tests,1));

% Test loop
for c = 1:num_tests
    
    close all
    
    % Set default parameters
    DefaultParameters;
    StructTemplates;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set test-specific parameters                                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Par.FLAG_AlgType = 0;
    Par.L = 1;
    Par.S = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set random stream
    s = RandStream('mt19937ar', 'seed', c);
    RandStream.setDefaultStream(s);
    
    % Generate target states
    [TrueTracks, InitStates] = GenerateStates();
    
    % Generate observations from target states
    [Observs, TrueTracks, detections] = GenerateObservations(TrueTracks);
    
%     % Plot states and observations
%     state_fig = PlotTrueTracks(TrueTracks);
%     obs_fig = PlotObs(Observs, detections);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Run tracking algorithm                                          %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if Par.FLAG_AlgType == 0
        [ Results ] = Track_MCMC(detections, Observs, InitStates );
    elseif Par.FLAG_AlgType == 1
        [ Results ] = Track_SISR(detections, Observs, InitStates );
    elseif Par.FLAG_AlgType == 2
        [ Results ] = Track_PDAF(detections, Observs, InitStates );
    elseif Par.FLAG_AlgType == 3
        [ Results ] = Track_JPDAF(detections, Observs, InitStates );
    elseif Par.FLAG_AlgType == 4
        [ Results ] = Track_MCMC_IS(detections, Observs, InitStates );
    elseif Par.FLAG_AlgType == 5
        [ Results ] = Track_UPDAF(detections, Observs, InitStates );
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot and analyse                                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1) || (Par.FLAG_AlgType == 4)
%         PlotParticles(Results{Par.T}.particles, state_fig);
%         [ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );
%     elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3) || (Par.FLAG_AlgType == 5)
%         PlotParticles(Results(Par.T), state_fig);
%     end
    
    if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1) || (Par.FLAG_AlgType == 4)
        [Stats.RMSE_MMSE(c), Stats.RMSE_MAP(c), Stats.prop_ass(c), Stats.num_lost(c)] = BasicParticleAnalysis(TrueTracks, Results);
        Stats.RMSE_MMSE(isnan(Stats.RMSE_MMSE))=0;
        total_RMSE = sqrt( sum(Stats.RMSE_MMSE.^2 .* (Par.NumTgts - Stats.num_lost))/sum(Par.NumTgts - Stats.num_lost) );
        total_prop_ass = sum( Stats.prop_ass .* (Par.NumTgts - Stats.num_lost) )/sum(Par.NumTgts - Stats.num_lost);
        total_lost = sum(Stats.num_lost)/(num_tests*Par.NumTgts);
    elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3)
        [Stats.RMSE(c), Stats.num_lost(c)] = BasicPointEstAnalysis(TrueTracks, Results);
        Stats.RMSE(isnan(Stats.RMSE))=0;
        total_RMSE = sqrt( sum(Stats.RMSE.^2 .* (Par.NumTgts - Stats.num_lost))/sum(Par.NumTgts - Stats.num_lost) );
        total_lost = sum(Stats.num_lost)/(num_tests*Par.NumTgts);
    end
    
end