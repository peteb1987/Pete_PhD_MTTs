% Run a single test with one algorithm.

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% Set default parameters
DefaultParameters;
StructTemplates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set test-specific parameters                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Par.FLAG_RB = true;
% Par.rand_seed = 2;
% Par.NumIt = 5000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set random stream
s = RandStream('mt19937ar', 'seed', Par.rand_seed);
RandStream.setDefaultStream(s);

% Generate target states
[TrueTracks, InitStates] = GenerateStates();

% Generate observations from target states
[Observs, detections] = GenerateObservations(TrueTracks);

% Plot states and observations
state_fig = PlotTrueTracks(TrueTracks);
obs_fig = PlotObs(Observs, detections);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run tracking algorithm                                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [ Results ] = Track_MCMC(detections, Observs, InitStates );
[ Results ] = Track_SISR(detections, Observs, InitStates );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot and analyse                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PlotParticles(Results{Par.T}.particles, state_fig);

[ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );