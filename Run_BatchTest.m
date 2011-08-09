% Run a batch of tests and analyse performance.

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% How many tests?
num_tests = 10;

Stats = repmat(struct('MMSE', [], 'LostTracks', []), num_tests, 1);

% Test loop
for c = 1:num_tests
    
    % Set default parameters
    DefaultParameters;
    StructTemplates;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set test-specific parameters                                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Par.FLAG_AlgType = 0;
    
    Par.T = 20;
    Par.NumTgts = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set random stream
    s = RandStream('mt19937ar', 'seed', c);
    RandStream.setDefaultStream(s);
    
    % Generate target states
    [TrueTracks, InitStates] = GenerateStates();
    
    % Generate observations from target states
    [Observs, TrueTracks, detections] = GenerateObservations(TrueTracks);
    
    % Plot states and observations
    state_fig = PlotTrueTracks(TrueTracks);
    obs_fig = PlotObs(Observs, detections);
    
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
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot and analyse                                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1)
        PlotParticles(Results{Par.T}.particles, state_fig);
        [ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );
    elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3)
        PlotParticles(Results(Par.T), state_fig);
    end
    
    if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1)
        [Stats(c).RMSE_MMSE, Stats(c).RMSE_MAP, Stats(c).prop_ass, Stats(c).num_lost] = BasicParticleAnalysis(TrueTracks, Results);
    elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3)
        [Stats(c).RMSE, Stats(c).num_lost] = BasicPointEstAnalysis(TrueTracks, Results);
    end
    
end