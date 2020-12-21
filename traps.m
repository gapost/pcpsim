function [xdot, V2, V2N] = traps(x,t,Xc0,V2_0,K1,K2)

V2 = V2_0 - (Xc0 - x); % traps V2
V2N = Xc0 - x;         %cluster of V2N


xdot = -K1.* V2 .* x + K2 .* V2N ; 

endfunction