% Define a set of empty structures to use as classes.

global Templates;
global Par;

% Observations
obs = struct('r', [], 'N', []);             % r is an array of observations (one per row).
                                            % N is the number of observations
Templates.Observs = repmat(obs, Par.T, 1);

% Track
Templates.Track = struct('num', [], 'birth', [], 'death', [], 'state', [], 'assoc', [], 'covar', [], 'smooth', []);
                % num is the number of frames for which the track lasts.
                % birth is the first frame in which the track appears
                % death is the first frame in which the track disappears
                % state is a cell array of state vectors, one for each frame
                % assoc is a vector of association indices, one for each frame
                % covar is a cell array of covariance matrices, one for each frame, when Gaussian approximations are being used (empty otherwise)

% TrackSet
Templates.TrackSet = struct('tracks', [], 'N', 0, 'members', [], 'origin', [], 'origin_time', []);
                % tracks is a cell array of Tracks
                % N is the number of tracks
                % members is the an arry of labels of the targets present, in the order they are in tracks
                % origin is the index of the particle from which the history is taken
                % origin_time is the processing frame from which the history was taken
                    % (these last two mainly for MCMC-PFs)

%%% Some useful lines of code for acting on these "classes"
% (Tracks{j}.birth <= t)&&(Tracks{j}.death > t)
% Tracks{j}.state{t-Tracks{j}.birth+1};