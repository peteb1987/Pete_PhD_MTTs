% Run a batch of tests and analyse performance.

% Clear the workspace (maintaining breakpoints)
clup

% Get test flag
test_flag = str2double(getenv('SGE_TASK_ID'))
% test_flag = 2

% How many tests?
num_tests = 11;

scaling = 0:0.1:1;

% Test loop
for c = 1:num_tests
    
    close all
    
    % Set default parameters
    DefaultParameters;
    StructTemplates;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set test-specific parameters                                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     switch mod(test_flag,8)
%         case 1
%             Par.FLAG_AlgType = 2;
%         case 2
%             Par.FLAG_AlgType = 3;
%         case 3
%             Par.FLAG_AlgType = 1;
%             Par.L = 1;
%         case 4
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%         case 5
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%             Par.FLAG_RB = true;
%         case 6
%             Par.FLAG_AlgType = 0;
%             Par.L = 1;
%             Par.S = 1;
%         case 7
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
%             Par.S = 5;
%         case 0
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
%             Par.S = 5;
%             Par.FLAG_RB = true;
%     end
%     
%     % Choose observation model
%     if test_flag > 8
%         Par.FLAG_ObsMod = 1;
%         Par.UnifPosDens = 1/(pi*Par.Xmax^2);
%         Par.ClutDens = (1/Par.Xmax)*(1/(2*pi));
%         Par.MinInitStateRadius = 0.25;
%         Par.MaxInitStateRadius = 0.35;
%         Par.BearingNoiseVar = 1E-4;
%         Par.RangeNoiseVar = 1;
%         Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
%     end
%     
%     Par.PDetect = 0.9;
%     Par.ExpClutObs = 200;
%     Par.FLAG_SetInitStates = true;
%     Par.InitStates = {[0; 200; 2; 0];
%         [0; 210; 2; 0];
%         [0; 220; 2; 0];
%         [0; 230; 2; 0];
%         [0; 240; 2; 0];};

    
%     % Choose algorithm
%     switch mod(test_flag,8)
%         case 1
%             Par.FLAG_AlgType = 2;
%         case 2
%             Par.FLAG_AlgType = 3;
%         case 3
%             Par.FLAG_AlgType = 1;
%             Par.L = 1;
%         case 4
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%         case 5
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%             Par.FLAG_RB = true;
%         case 6
%             Par.FLAG_AlgType = 0;
%             Par.L = 1;
%             Par.S = 1;
%         case 7
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
%             Par.S = 5;
%         case 0
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
%             Par.S = 5;
%             Par.FLAG_RB = true;
%     end
%     
%     % Choose observation model
%     if test_flag > 32
%         Par.FLAG_ObsMod = 1;
%         Par.UnifPosDens = 1/(pi*Par.Xmax^2);
%         Par.ClutDens = (1/Par.Xmax)*(1/(2*pi));
%         Par.MinInitStateRadius = 0.25;
%         Par.MaxInitStateRadius = 0.35;
%         Par.BearingNoiseVar = 1E-4;
%         Par.RangeNoiseVar = 1;
%         Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
%     end
%     
%     % Choose scenario
%     switch ceil(test_flag/8)
%         case {1,5}
%             
%         case {2,6}
%             Par.PDetect = 0.9;
%             Par.ExpClutObs = 200;
%             
%         case {3,7}
%             Par.PDetect = 0.9;
%             Par.ExpClutObs = 200;
%             Par.ProcNoiseVar = 10;
%             Par.Q = Par.ProcNoiseVar * [P^3/3 0 P^2/2 0; 0 P^3/3 0 P^2/2; P^2/2 0 P 0; 0 P^2/2 0 P];
%             if Par.FLAG_ObsMod == 0
%                 Par.ObsNoiseVar = 10;
%                 Par.R = Par.ObsNoiseVar * eye(2);
%             elseif Par.FLAG_ObsMod == 1
%                 Par.BearingNoiseVar = 1E-3;
%                 Par.RangeNoiseVar = 10;
%                 Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
%             end
%             
%         case {4,8}
%             Par.FLAG_SetInitStates = true;
%             Par.InitStates = {[0; 200; 2; 0];
%                 [0; 210; 2; 0];
%                 [0; 220; 2; 0];
%                 [0; 230; 2; 0];
%                 [0; 240; 2; 0];};
% 
%               
%     end

    % Choose algorithm
    Par.FLAG_AlgType = 4;
    
    % Set acceptance scaling
    Par.HistoryAcceptScaling = scaling(test_flag);
    
%     % Choose observation model
%     if test_flag > 4
%         Par.FLAG_ObsMod = 1;
%         Par.UnifPosDens = 1/(pi*Par.Xmax^2);
%         Par.ClutDens = (1/Par.Xmax)*(1/(2*pi));
%         Par.MinInitStateRadius = 0.25;
%         Par.MaxInitStateRadius = 0.35;
%         Par.BearingNoiseVar = 1E-4;
%         Par.RangeNoiseVar = 1;
%         Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
%     end
%     
%     % Choose scenario
%     switch test_flag
%         case {1,5}
%             
%         case {2,6}
%             Par.PDetect = 0.9;
%             Par.ExpClutObs = 200;
%             
%         case {3,7}
%             Par.PDetect = 0.9;
%             Par.ExpClutObs = 200;
%             Par.ProcNoiseVar = 10;
%             Par.Q = Par.ProcNoiseVar * [P^3/3 0 P^2/2 0; 0 P^3/3 0 P^2/2; P^2/2 0 P 0; 0 P^2/2 0 P];
%             if Par.FLAG_ObsMod == 0
%                 Par.ObsNoiseVar = 10;
%                 Par.R = Par.ObsNoiseVar * eye(2);
%             elseif Par.FLAG_ObsMod == 1
%                 Par.BearingNoiseVar = 1E-3;
%                 Par.RangeNoiseVar = 10;
%                 Par.R = [Par.BearingNoiseVar 0; 0 Par.RangeNoiseVar];
%             end
%             
%         case {4,8}
%             Par.FLAG_SetInitStates = true;
%             Par.InitStates = {[0; 200; 2; 0];
%                 [0; 210; 2; 0];
%                 [0; 220; 2; 0];
%                 [0; 230; 2; 0];
%                 [0; 240; 2; 0];};
% 
%               
%     end
    
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
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot and analyse                                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1)
%         PlotParticles(Results{Par.T}.particles, state_fig);
%         [ assoc ] = RetrieveAssocs( Par.T, Results{Par.T}.particles );
%     elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3)
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
        total_prop_ass = 0;
    end
    
end

save(['test' num2str(test_flag) '.mat'], 'Par', 'total_RMSE', 'total_lost', 'total_prop_ass', 'Stats');