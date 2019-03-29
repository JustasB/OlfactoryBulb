%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function set the cluster paramters
%
% Licurgo de Almeida
% 11/20/2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,minworker,maxworker] = SetCluster(name)
if strcmp(name,'GPU')
    time = 300;
    minworker = 8;
    maxworker = 8;
elseif strcmp(name,'Quick')
    time = 10;
    minworker = 4;
    maxworker = 4;
elseif strcmp(name,'Default')
    time = 300;
    minworker = 8;
    maxworker = 51;
end