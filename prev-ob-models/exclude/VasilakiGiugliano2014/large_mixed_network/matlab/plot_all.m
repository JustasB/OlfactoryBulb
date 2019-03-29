%
%
%
clear all;
clc;

figure(1); set(gcf, 'DoubleBuffer', 'on');

addpath matlab;
%addpath Gen_and_Analize_Connect;
P = dir('W_*.x');
N = length(P);
load connectivity.mat;

fname = P(end).name
%W     = load(fname);

H = mysubplot(3,2,0.1,0.1,0.05);

axes(H(1));
%weight_distribution(W);

axes(H(2));
plot_rate;

axes(H(3));
plot_rate_histogram;

%axes(H(4));
%plot_raster;

% axes(H(5));
% [out nil p] = statistics_motifs_pairs_normal(abs(C), W, 5., 2./3.);
% %[out nil p] = statistics_motifs_pairs_normal(abs(C), abs(C), 1, 0.001);
% 
% axes(H(6));
% [out nil p] = statistics_motifs_pairs_facildepress(C, abs(C), W, 5., 2./3.);

set(H, 'FontSize', 12);

