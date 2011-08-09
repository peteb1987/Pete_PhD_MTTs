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
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Set random stream
    s = RandStream('mt19937ar', 'seed', c);
    RandStream.setDefaultStream(s);
    
    % Generate target states
    [TrueTracks, InitStates] = GenerateStates();
    
    % Generate observations from target states
    [Observs, detections] = GenerateObservations(TrueTracks);
    
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
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot and analyse                                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Stats(c).MSE, Stats(c).LostTracks] = BasicParticleAnalysis(TrueTracks, Results);
    
end