% This adds the current directory ("ModelDB") and all subdirectories to the
% path. Essential for much of the included code to run properly.

[pathstr, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(pathstr))
clear pathstr