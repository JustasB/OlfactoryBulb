%
%
%
%

clear all;
close all;
clc;

N = 500;        % Size of the simulated network whose connectivity must be
                % synthetically generated...

%--------------------------------------------------------------------------
% Test 1 - random connectivity
p = 0.4;        % Marginal probability of a synaptic connection C(i,j).

P1 = 2*p*(1-p); % Resulting probability of unidirectional motifs.
P2 = p*p;       % Resulting probability of bidirectional motifs.
C1 = generate_connectivity(N, P1, P2); % Generate connectivity

% Let's evaluate the resulting connectivity matrix, analyzing its motifs
% according to the actual data and against the null hypothesis of random
% connectivity
Amax = 1;       % Parameters set to make the script working with C = {0,1}
h    = 0.001;   % Parameters set to make the script working with C = {0,1}
[out1 nil1 p_estimate1] = statistics_motifs_pairs_normal(C1, C1, Amax, h);

disp(sprintf('Connectivity is by construction completely random.'));
disp(sprintf('i.e. NULL hypothesis cannot be rejected.\n'));
disp(sprintf('p = %f, estimated as %f', p, p_estimate1));
disp(sprintf('Motifs probs: %f %f %f', out1(1), out1(2), out1(3)));
disp(sprintf('Null hypoth : %f %f %f\n\n', nil1(1), nil1(2), nil1(3)));

figure(1);
disp(sprintf('Press enter for the next test...\n\n'));
pause;

%--------------------------------------------------------------------------
% Test 2 - non-random connectivity

P1 = 0.4;
P2 = 0.4;
C2 = generate_connectivity(N, P1, P2); % Generate connectivity

% Let's evaluate the resulting connectivity matrix, analyzing its motifs
% according to the actual data and against the null hypothesis of random
% connectivity
Amax = 1;       % Parameters set to make the script working with C = {0,1}
h    = 0.001;   % Parameters set to make the script working with C = {0,1}
[out2 nil2 p_estimate2] = statistics_motifs_pairs_normal(C2, C2, Amax, h);

disp(sprintf('Connectivity is by construction NON-random.'));
disp(sprintf('i.e. NULL hypothesis is rejected.\n'));
disp(sprintf('Unidirectional and bidirectional motifs should have same prob.\n'));
disp(sprintf('Motifs probs: %f %f %f', out2(1), out2(2), out2(3)));
disp(sprintf('Null hypoth : %f %f %f\n\n', nil2(1), nil2(2), nil2(3)));

figure(1);
disp(sprintf('Press enter for the next test...\n\n'));
pause;

%--------------------------------------------------------------------------


% Test 3 - random connectivity and short-term dynamics
Pff= 0.25;
Pdd= 0.25;
Pfd= 0.25;
%Pfd = Pfd
Pf = 0.5;
Pd = 0.5;

Dyn = generate_dynamics(C1, Pff, Pdd, Pfd, Pf, Pd);
[out nil p] = statistics_motifs_pairs_facildepress(Dyn, C1, C1, Amax, h);

disp(sprintf('Connectivity is by construction random.'));
disp(sprintf('i.e. NULL hypothesis cannot be rejected.\n'));
disp(sprintf('Motifs have equal probabilities of occurrence\n'));
disp(sprintf('Motifs probs: %f %f %f %f %f %f', out(1), out(2), out(3), out(4), out(5), out(6)));
disp(sprintf('Null hypoth : %f %f %f %f %f %f\n\n', nil(1), nil(2), nil(3), nil(4), nil(5), nil(6)));

figure(1);
disp(sprintf('Press enter for the next test...\n\n'));
pause;

%--------------------------------------------------------------------------


% Test 4 - random connectivity and short-term dynamics
Pff= 0.1;
Pdd= 0.8;
Pfd= 0.05;
%Pfd = Pfd
Pf = 0.75;
Pd = 0.25;

Dyn = generate_dynamics(C1, Pff, Pdd, Pfd, Pf, Pd);
[out nil p] = statistics_motifs_pairs_facildepress(Dyn, C1, C1, Amax, h);

disp(sprintf('Connectivity is by construction random but synaptic dynamics is NOT.'));
disp(sprintf('i.e. NULL hypothesis can be rejected.\n'));
disp(sprintf('Motifs have different probabilities of occurrence\n'));
disp(sprintf('Motifs probs: %f %f %f %f %f %f', out(1), out(2), out(3), out(4), out(5), out(6)));
disp(sprintf('Null hypoth : %f %f %f %f %f %f\n\n', nil(1), nil(2), nil(3), nil(4), nil(5), nil(6)));

figure(1);