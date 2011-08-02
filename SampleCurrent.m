function [ppsl_prob, assoc, state] = SampleCurrent(j, t, L, Set, Observs, no_samp)
%SAMPLECURRENT Sample the current associations and states

[assoc, assoc_prob] = SampleAssociations(j, t, L, Set, Observs, no_samp);
Set.tracks(j).assoc(t-L+1 -Set.tracks(j).birth+1:t -Set.tracks(j).birth+1) = assoc;
[state, state_prob] = SampleStates(j, t, L, Set, Observs, no_samp);
ppsl_prob = assoc_prob + state_prob;

end

