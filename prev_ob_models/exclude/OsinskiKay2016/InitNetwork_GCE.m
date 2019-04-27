
function [dt,tsim,numtp,nmit,ngradist,ngraprox,sampf,timevec] = InitNetwork_GCE(input_file)

% Boleslaw Osinski (2015)

[crap, nums] = textread(input_file,'%s %s',40);
dt = str2num(nums{2});
tsim = str2num(nums{3});
numtp = length(1:round(tsim/dt));
nmit = str2num(nums{11});
ngradist = str2num(nums{12});
ngraprox = str2num(nums{13});

sampf = 1/(dt*1e-3);

timevec = dt*(1:numtp);

end