% Run a single test with one algorithm.

% Clear the workspace (maintaining breakpoints)
clup
dbstop if error

% Set default parameters
DefaultParameters;
StructTemplates;

% Set random stream
s = RandStream('mt19937ar', 'seed', Par.rand_seed);
RandStream.setDefaultStream(s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set test-specific parameters                                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate target states
[TrueTracks, InitStates] = GenerateStates();

% Generate observations from target states
[Observs, detections] = GenerateObservations(TrueTracks);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run tracking algorithm                                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ Results ] = Track_MCMC(detections, Observs, InitStates );







% Plot states and observations
fig_states = PlotTrueState(TrueState);
PlotObs(Observs, detections);

% Run tracker
[ Chains ] = MultiTargetTrack(detections, Observs, {TargSpec(:).state} );

% Plot final estimates
PlotTracks(Chains{Par.T}, fig);

% Analyse associations
[ass, count, present] = AnalyseAss( detections, Chains{Par.T}, Par.T);