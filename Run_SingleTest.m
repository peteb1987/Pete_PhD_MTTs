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

% Par.rand_seed = 2;
Par.T = 20;
Par.NumTgts = 1;
% Par.FLAG_SetInitStates = true;
% Par.InitStates = {[0; 200; 2; 0];
%                   [0; 210; 2; 0];
%                   [0; 220; 2; 0];
%                   [0; 230; 2; 0];
%                   [0; 240; 2; 0];};
% Par.ProcNoiseVar = 10;
% Par.Q = Par.ProcNoiseVar * [P^3/3 0 P^2/2 0; 0 P^3/3 0 P^2/2; P^2/2 0 P 0; 0 P^2/2 0 P];
% Par.ObsNoiseVar = 10;
% Par.R = Par.ObsNoiseVar * eye(2);
% Par.ExpClutObs = 200;
% Par.PDetect = 0.9;

% Par.BearingNoiseVar = 1E-3;
% Par.RangeNoiseVar = 10;
% Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];

% Par.rand_seed = 2;

Par.FLAG_AlgType = 4;
% Par.FLAG_DynMod = 1;
% Par.FLAG_ObsMod = 1;
% Par.L = 5;
% Par.S = 5;

% Par.FLAG_RB = true;
Par.NumIt = 500;

% Par.NumTgts = 1;
% Par.T = 10;

% Par.PDetect = 0.9;
% Par.ExpClutObs = 200;
% Par.ProcNoiseVar = 10;
% Par.Q = Par.ProcNoiseVar * [P^3/3 0 P^2/2 0; 0 P^3/3 0 P^2/2; P^2/2 0 P 0; 0 P^2/2 0 P];
% if Par.FLAG_ObsMod == 0
%     Par.ObsNoiseVar = 10;
%     Par.R = Par.ObsNoiseVar * eye(2);
% elseif Par.FLAG_ObsMod == 1
%     Par.BearingNoiseVar = 1E-3;
%     Par.RangeNoiseVar = 10;
%     Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
% end

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot and analyse                                                    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1) || (Par.FLAG_AlgType == 4)
    PlotParticles(Results{Par.T}.particles, state_fig);
    [ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );
elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3)
    PlotParticles(Results(Par.T), state_fig);
end