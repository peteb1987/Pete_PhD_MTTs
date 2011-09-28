function [ state ] = IntrinsicDynamicInverse( next_state, acc )
%INTRINSICDYNAMICINVERSE Evaluate inverse state transition


global Par;
P = Par.P;

aT = acc(1);
aP = acc(2);
ax1 = acc(3);
ax2 = acc(4);

state = zeros(4,1);

new_phi = next_state(3);
new_sdot = next_state(4);
new_x1 = next_state(1);
new_x2 = next_state(2);

SF1 = 4*aT^2 + aP^2;

sdot = new_sdot - aT*P;
if sdot<Par.MinSpeed
    sdot=Par.MinSpeed;
    aT = (new_sdot-sdot)/P;
end
if aT~=0
    phi = new_phi - (aP/aT)*log(new_sdot/sdot);
else
    phi = new_phi - (aP*P)/sdot;
end

if (aT~=0)&&(aP~=0)
    state(1) = new_x1 - ((new_sdot^2/SF1)*( aP*sin(new_phi)+2*aT*cos(new_phi)) - (sdot^2/SF1)*( aP*sin(phi)+2*aT*cos(phi))) - ax1;
    state(2) = new_x2 - ((new_sdot^2/SF1)*(-aP*cos(new_phi)+2*aT*sin(new_phi)) - (sdot^2/SF1)*(-aP*cos(phi)+2*aT*sin(phi))) - ax2;
elseif (aT==0)&&(aP~=0)
    state(1) = new_x1 - (new_sdot^2/aP)*( sin(new_phi) - sin(phi) ) - ax1;
    state(2) = new_x2 - (new_sdot^2/aP)*( -cos(new_phi) + cos(phi) ) - ax2;
elseif (aT~=0)&&(aP==0)
    state(1) = new_x1 - 0.5*Par.P*cos(phi)*new_sdot - ax1;
    state(2) = new_x2 - 0.5*Par.P*sin(phi)*new_sdot - ax2;
else
    state(1) = new_x1 - ax1 - ( sdot*P*cos(phi) );
    state(2) = new_x2 - ax2 - ( sdot*P*sin(phi) );
end
state(3) = phi;
state(4) = sdot;

% assert(~any(isnan(state)), 'NaN state!');


end

