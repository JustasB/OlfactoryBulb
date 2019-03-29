%
% Launcher application for the Large Networks C-code simulations
% (first connectivity is generated and then the simulation launched)
%

clear all;      % clear the workspace
close all;      % close any open figure or file
clc;            % clear the screen
addpath matlab; % add path to relevant routines in the ./matlab directory


Ncell = 1000;    % Number of excitatory cells constituting the network 


fixseed   = 0;                          % currently unused
seed      = 3532765;                    % currently unused
pars      = [2 350 150 1000 160 3];     % currently unused


%--------------------------------------------------------------------------
% Stat. description of  structural, hard-wired, anatomical connectivity
%
Pdd       = .8;
Pff       = .8;
Pfd       = .8;
Pdf       = .8;
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    StructuredMatrix(Ncell, Pdd, Pff, Pfd, Pdf);
%--------------------------------------------------------------------------
%save connectivity information
%-----------------------------
% % % % % % % % % % % % % % MM    = length(find(C(:)~=0));  % Number of non-zero elements of C(i,j)
% % % % % % % % % % % % % % Ctemp = C';
% % % % % % % % % % % % % % Ctemp = Ctemp(:);               % A single vector, encoding 'C' row-wise
% % % % % % % % % % % % % % tmp   = find(Ctemp ~= 0);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % fname = 'connectivity.dat';
% % % % % % % % % % % % % % fp    = fopen(fname, 'wb');
% % % % % % % % % % % % % % fwrite(fp, Ncell, 'double');  
% % % % % % % % % % % % % % fwrite(fp, MM, 'double');  
% % % % % % % % % % % % % % %for i=1:Ncell, fwrite(fp, int8(C(i,:)), 'double'); end % it was 'short'
% % % % % % % % % % % % % % for kk=1:length(tmp)
% % % % % % % % % % % % % %     fwrite(fp, tmp(kk) * sign(Ctemp(tmp(kk))), 'double');
% % % % % % % % % % % % % % end
% % % % % % % % % % % % % % fclose(fp);
% % % % % % % % % % % % % % 
% % % % % % % % % % % % % % save('connectivity.mat', 'C');
% % % % % % % % % % % % % % disp('Connectivity matrix generation and storage on disk completed!');
%--------------------------------------------------------------------------
% save some human-redable information on file
% % % % % % % % % % fname = 'info.txt';
% % % % % % % % % % fp    = fopen(fname, 'w');
% % % % % % % % % % fprintf(fp,'Simulation input details and parameters\n\n');
% % % % % % % % % % fprintf(fp,'Ncell = %d\n', Ncell);
% % % % % % % % % % fprintf(fp,'seed (fix = %d)  = %d\n', fixseed, seed);
% % % % % % % % % % fprintf(fp,'---- Probabilities ----\n');
% % % % % % % % % % fprintf(fp,'[Prob_ee Prob_loop]\t\t\t= [%.3f\t%.3f]\n', Prob_ee, Prob_loop);
% % % % % % % % % % fprintf(fp,'[eeff eedd eefd eedf]\t\t= [%.3f\t%.3f\t%.3f\t%.3f]\n', Prob_ff , Prob_dd, Prob_fd, Prob_fd);
% % % % % % % % % % fprintf(fp,'[exf exd]\t\t= [%.3f\t%.3f]\n', Prob_f, Prob_d);
% % % % % % % % % % fclose(fp);
% % % % % % % % % % %-------------------
% % % % % % % % % % %all_pars = [pars, Ncell, Prob_ee, Prob_loop, Prob_ff, Prob_dd, Prob_fd, Prob_f, Prob_d, fixseed, seed];
% % % % % % % % % % all_pars = [Ncell, Prob_ee, Prob_loop, Prob_ff, Prob_dd, Prob_fd, Prob_f, Prob_d, fixseed, seed];
% % % % % % % % % % fname = 'pars.dat';
% % % % % % % % % % fp    = fopen(fname, 'wb');
% % % % % % % % % % fwrite(fp, pars, 'double');  
% % % % % % % % % % fclose(fp);
% % % % % % % % % % save('pars.mat', 'all_pars');
%--------------------------------------------------------------------------


% Usage: IFnet Tsim [sec]  Io [pA] GA [pA] (was 200.)
% cmd = sprintf('./IFnet 100 200 12');
% disp(sprintf('Starting the simulation [%s]...', cmd));
% [status,result] = system(cmd)


