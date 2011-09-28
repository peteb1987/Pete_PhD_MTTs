% Run a single test with one algorithm.

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error
% dbstop if warning

% Set default parameters
DefaultParameters;
StructTemplates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set test-specific parameters                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.rand_seed = 9;
Par.FLAG_AlgType = 0;
% Par.FLAG_RB = true;
Par.L = 5;
Par.S = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random stream
s = RandStream('mt19937ar', 'seed', Par.rand_seed);
RandStream.setDefaultStream(s);

% Generate target states
[TrueTracks, InitStates] = GenerateStates();

% Generate observations from target states
[Observs, TrueTracks, detections] = GenerateObservations(TrueTracks);

% Plot states and observations
state_fig = PlotTrueTracks(TrueTracks);
obs_fig = PlotObs(Observs, detections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run tracking algorithm                                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot and analyse                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1) || (Par.FLAG_AlgType == 4)
    PlotParticles(Results{Par.T}.particles, state_fig);
    [ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );
elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3) || (Par.FLAG_AlgType == 5)
    PlotParticles(Results(Par.T), state_fig);
end