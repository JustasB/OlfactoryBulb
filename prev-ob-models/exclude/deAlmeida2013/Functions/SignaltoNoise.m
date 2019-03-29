function S = SignaltoNoise(C1,C2)
% function S = SignaltoNoise(C1,C2) calculates the signal to noise ratio in
% the network. Where:
%
if isobject(C1)
    a1 = sum(sum(C1.O));
    a2 = sum(sum(C2.O));
else
    a1 = sum(C1);
    a2 = sum(C2);
end
%S = a1 / (a1 + a2);
S = 1 - (a2 / a1);
if S > 1
    S = 1;
elseif S < 0
    S = 0;
end
end