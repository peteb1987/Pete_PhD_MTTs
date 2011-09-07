function [ A, Q ] = IntrinsicDynamicLinearise( state )
%INTRINSICDYNAMICLINEARISE Calculate linearised A and Q matrices for a
% given state for the intrinsic coordinate dynamic model

global Par;

P = Par.P;

phi = state(3); sdot = state(4);
A = [1 0 -P*sdot*sin(phi) P*cos(phi);
    0 1   P*sdot*cos(phi) P*sin(phi);
    0 0   1               0         ;
    0 0   0               1          ];

B = [0.5*P^2*cos(phi) -0.5*P^2*sin(phi) 1 0;
    0.5*P^2*sin(phi)   0.5*P^2*cos(phi) 0 1;
    0                  P/sdot           0 0;
    P                  0                0 0 ];

Q = B*Par.Q_pre*B';

end

