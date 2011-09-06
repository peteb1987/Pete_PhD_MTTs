function [ state ] = IntrinsicDynamicEvaluate( prev_state, aT, aP )
%INTRINSICDYNAMICEVALUATE Calculate the next state value given the random
% accelerations.

global Par;
P = Par.P;

state = zeros(4,1);

phi = prev_state(3);
sdot = prev_state(4);
x = prev_state(1);
y = prev_state(2);

SF1 = 4*aT^2 + aP^2;
SF2 = (aT*P+sdot)^2;

new_sdot = max(sdot + aT*P, 0.1);
new_phi = phi + (aP/aT)*log(new_sdot/sdot);

state(1) = x + (SF2/SF1)*( aP*sin(new_phi)+2*aT*cos(new_phi)) - (sdot^2/SF1)*( aP*sin(phi)+2*aT*cos(phi));
state(2) = y + (SF2/SF1)*(-aP*cos(new_phi)+2*aT*sin(new_phi)) - (sdot^2/SF1)*(-aP*cos(phi)+2*aT*sin(phi));
state(3) = new_phi;
state(4) = new_sdot;

end

