function [ track ] = ProjectTrack( t, track )
%PROJECTTRACK Extend and project track forward using ML estimate

global Par;

prev_state = track.state{t-1-track.birth+1};

if Par.FLAG_DynMod == 0
    state = Par.A * prev_state;
    covar = Par.A * track.covar{t-1-track.birth+1} * Par.A + Par.Q;
elseif Par.FLAG_DynMod == 1
    [A, Q] = IntrinsicDynamicLinearise(prev_state);
    state = IntrinsicDynamicPredict(prev_state);
    covar = A * track.covar{t-1-track.birth+1} * A + Q;
end

track.death = track.death + 1;
track.num = track.num + 1;
track.state = [track.state; {state}];
track.smooth = [track.smooth; {state}];
track.covar = [track.covar; {covar}];
track.assoc = [track.assoc; 0];

end

