function [ state ] = IntrinsicDynamicEvaluate( prev_state, acc )
%INTRINSICDYNAMICEVALUATE Calculate the next state value given the random
% accelerations.

global Par;
P = Par.P;

aT = acc(1);
aP = acc(2);
ax1 = acc(3);
ax2 = acc(4);

state = zeros(4,1);

phi = prev_state(3);
sdot = prev_state(4);
x1 = prev_state(1);
x2 = prev_state(2);

SF1 = 4*aT^2 + aP^2;
SF2 = (aT*P+sdot)^2;

new_sdot = sdot + aT*P;
if new_sdot<Par.MinSpeed
    new_sdot=Par.MinSpeed;
    aT = (new_sdot-sdot)/P;
end
if aT~=0
    new_phi = phi + (aP/aT)*log(new_sdot/sdot);
else
    new_phi = phi + (aP*P)/sdot;
end

state(1) = x1 + (SF2/SF1)*( aP*sin(new_phi)+2*aT*cos(new_phi)) - (sdot^2/SF1)*( aP*sin(phi)+2*aT*cos(phi)) + ax1;
state(2) = x2 + (SF2/SF1)*(-aP*cos(new_phi)+2*aT*sin(new_phi)) - (sdot^2/SF1)*(-aP*cos(phi)+2*aT*sin(phi)) + ax2;
state(3) = new_phi;
state(4) = new_sdot;

assert(~any(isnan(state)), 'NaN state!');

end

