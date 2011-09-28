function [ state ] = IntrinsicDynamicEvaluate( prev_state, acc )
%INTRINSICDYNAMICEVALUATE Calculate the next state value given the random
% accelerations.

global Par;
min_speed = Par.MinSpeed;
P = Par.P;

Np = size(acc, 2);

aT = acc(1,:);
aP = acc(2,:);
ax1 = acc(3,:);
ax2 = acc(4,:);

state = zeros(4,Np);

phi = prev_state(3,:);
sdot = prev_state(4,:);
x1 = prev_state(1,:);
x2 = prev_state(2,:);

new_sdot = sdot + aT*P;
new_sdot(new_sdot<min_speed) = min_speed;
aT = (new_sdot-sdot)/P;

SF1 = 4*aT.^2 + aP.^2;
SF2 = new_sdot.^2;

if aT~=0
    new_phi = phi + (aP./aT).*log(new_sdot./sdot);
else
    new_phi = phi + (aP.*P)./sdot;
end


for ii = 1:Np
    if (aT(ii)~=0)&&(aP(ii)~=0)
        state(1,ii) = x1(ii) + (SF2(ii)./SF1(ii)).*( aP(ii).*sin(new_phi(ii))+2*aT(ii).*cos(new_phi(ii))) - (sdot(ii)^2./SF1(ii))*( aP(ii).*sin(phi(ii))+2*aT(ii).*cos(phi(ii))) + ax1(ii);
        state(2,ii) = x2(ii) + (SF2(ii)./SF1(ii)).*(-aP(ii).*cos(new_phi(ii))+2*aT(ii).*sin(new_phi(ii))) - (sdot(ii)^2./SF1(ii))*(-aP(ii).*cos(phi(ii))+2*aT(ii).*sin(phi(ii))) + ax2(ii);
    elseif (aT(ii)==0)&&(aP(ii)~=0)
        state(1,ii) = x1(ii) + (SF2(ii)./aP(ii)).*( sin(new_phi(ii)) - sin(phi(ii)) ) + ax1(ii);
        state(2,ii) = x2(ii) + (SF2(ii)./aP(ii)).*( -cos(new_phi(ii)) + cos(phi(ii)) ) + ax2(ii);
    elseif (aT(ii)~=0)&&(aP(ii)==0)
        state(1,ii) = x1(ii) + 0.5*P*cos(phi(ii)).*new_sdot(ii) + ax1(ii);
        state(2,ii) = x2(ii) + 0.5*P*sin(phi(ii)).*new_sdot(ii) + ax2(ii);
    else
        state(1,ii) = x1(ii) + ax1(ii) + ( sdot(ii)*P.*cos(phi(ii)) );
        state(2,ii) = x2(ii) + ax2(ii) + ( sdot(ii)*P.*sin(phi(ii)) );
    end
end
state(3,:) = new_phi;
state(4,:) = new_sdot;

% assert(~any(isnan(state)), 'NaN state!');
% assert(all(isreal(state)), 'Complex state!');

end

