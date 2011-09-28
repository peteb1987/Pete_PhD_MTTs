% Run a batch of tests and analyse performance.

% Clear the workspace (maintaining breakpoints)
clup

% Get test flag
test_flag = str2double(getenv('SGE_TASK_ID'))

% How many tests?
num_tests = 20;

% Test loop
for c = 1:num_tests
    
    close all
    
    % Set default parameters
    DefaultParameters;
    StructTemplates;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set test-specific parameters                                    %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    title = 'MCMCIS_';
    
    curr = mod(test_flag-1,6)+1;
    hist = rem(test_flag-1,6)+1;
    scaling = 0:0.2:1;
    Par.CurrentAcceptScaling = scaling(curr);
    Par.HistoryAcceptScaling = scaling(hist);
    Par.FLAG_AlgType = 4;
    
%     % Choose algorithm
%     switch test_flag
%         case 1
%             Par.FLAG_AlgType = 2;
%         case 2
%             Par.FLAG_AlgType = 5;
%         case 3
%             Par.FLAG_AlgType = 1;
%             Par.L = 1;
%         case 4
%             Par.FLAG_AlgType = 1;
%             Par.L = 3;
%         case 5
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%         case 6
%             Par.FLAG_AlgType = 0;
%             Par.L = 1;
% %             Par.S = 1;
%         case 7
%             Par.FLAG_AlgType = 0;
%             Par.L = 3;
% %             Par.S = 3;
%         case 8
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
% %             Par.S = 5;
%         case 9
%             Par.FLAG_AlgType = 1;
%             Par.L = 1;
%             Par.FLAG_RB = true;
%         case 10
%             Par.FLAG_AlgType = 1;
%             Par.L = 3;
%             Par.FLAG_RB = true;
%         case 11
%             Par.FLAG_AlgType = 1;
%             Par.L = 5;
%             Par.FLAG_RB = true;
%         case 12
%             Par.FLAG_AlgType = 0;
%             Par.L = 1;
% %             Par.S = 1;
%             Par.FLAG_RB = true;
%         case 13
%             Par.FLAG_AlgType = 0;
%             Par.L = 3;
% %             Par.S = 3;
%             Par.FLAG_RB = true;
%         case 14
%             Par.FLAG_AlgType = 0;
%             Par.L = 5;
% %             Par.S = 5;
%             Par.FLAG_RB = true;
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
    elseif Par.FLAG_AlgType == 5
        [ Results ] = Track_UPDAF(detections, Observs, InitStates );
    end
    
%     save([title '_Results_test' num2str(test_flag) '_run' num2str(c) '.mat'], 'Results', 'Observs', 'TrueTracks', 'detections');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot and analyse                                                %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Par.FLAG_AlgType == 0) || (Par.FLAG_AlgType == 1) || (Par.FLAG_AlgType == 4)
        [Stats.RMSE_MMSE(c), Stats.RMSE_MAP(c), Stats.prop_ass(c), Stats.num_lost(c)] = BasicParticleAnalysis(TrueTracks, Results);
        Stats.RMSE_MMSE(isnan(Stats.RMSE_MMSE))=0;
        total_RMSE = sqrt( sum(Stats.RMSE_MMSE.^2 .* (Par.NumTgts - Stats.num_lost))/sum(Par.NumTgts - Stats.num_lost) );
        total_prop_ass = sum( Stats.prop_ass .* (Par.NumTgts - Stats.num_lost) )/sum(Par.NumTgts - Stats.num_lost);
        total_lost = sum(Stats.num_lost)/(num_tests*Par.NumTgts);
    elseif (Par.FLAG_AlgType == 2) || (Par.FLAG_AlgType == 3) || (Par.FLAG_AlgType == 5)
        [Stats.RMSE(c), Stats.num_lost(c)] = BasicPointEstAnalysis(TrueTracks, Results);
        Stats.RMSE(isnan(Stats.RMSE))=0;
        total_RMSE = sqrt( sum(Stats.RMSE.^2 .* (Par.NumTgts - Stats.num_lost))/sum(Par.NumTgts - Stats.num_lost) );
        total_lost = sum(Stats.num_lost)/(num_tests*Par.NumTgts);
        total_prop_ass = 0;
    end
    
end

save([title '_test' num2str(test_flag) '.mat'], 'Par', 'total_RMSE', 'total_lost', 'total_prop_ass', 'Stats');