function [ state ] = IntrinsicDynamicPredict( prev_state )
%INTRINSICDYNAMICPREDICT Project track forward 1 step

global Par;

state = zeros(4,1);
state(1) = prev_state(1) + prev_state(4)*Par.P*cos(prev_state(3));
state(2) = prev_state(2) + prev_state(4)*Par.P*sin(prev_state(3));
state(3) = prev_state(3);
state(4) = prev_state(4);

end

