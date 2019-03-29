function C = generate_dynamics(C, Pff, Pdd, Pfd, Pf, Pd)
%
% This function modifies the elements of the input binary connectivity 
% matrix [C], whose generic as it follows.
%
% Convention: C(i,j) == 1, iff facilitaiting connection
%             C(i,j) == -1, iff depressing connection
%
% if they are '1' they can be turned into '1' or '-1' with the specified set
% of probabilities...
%
% For those cases in which C(i,j) == 1 and C(j,i) == 1:
%
% Pff --> probability that C(i,j) will be set to '1' and that C(j,i).. '1'
% Pdd --> probability that C(i,j) will be set to '-1' and that C(j,i).. '-1'
% Pfd --> probability that C(i,j) will be set to '1' and that C(j,i).. '-1'
% Pdf --> probability that C(i,j) will be set to '-1' and that C(j,i).. '1'
%
% For the other cases
%
% Pf --> probability that C(i,j) will be set to '1' while that C(j,i) was '0'
% Pd --> probability that C(i,j) will be set to '-1' while that C(j,i) was '0'
%
%
%

eventlist_loop     = [1 2 3 4];
probabilities_loop = [Pff Pdd Pfd Pfd];

eventlist     = [5 6];
probabilities = [Pf Pd];

N = size(C,1);

for i=1:N,
  for j=(i+1):N,
      
    event = 0;
    
    if (C(i,j)  && C(j,i)),  event = return_event(eventlist_loop, probabilities_loop); end
    if (~C(i,j) && C(j,i)),  event = return_event(eventlist, probabilities); end
    if (C(i,j)  && ~C(j,i)), event = return_event(eventlist, probabilities); end
       
    switch event
        case 1,
            C(i,j) = 1;
            C(j,i) = 1;
        case 2,
            C(i,j) = -1;
            C(j,i) = -1;            
        case 3,
            C(i,j) = 1;
            C(j,i) = -1;            
        case 4,
            C(i,j) = -1;
            C(j,i) = 1;            
        case 5,
            C(i,j) = C(i,j) * 1;
            C(j,i) = C(j,i) * 1;
        case 6,
            C(i,j) = C(i,j) * -1;
            C(j,i) = C(j,i) * -1;
    end
    
    end
end


end







function event = return_event(eventlist, probabilities)
%
%
% e.g. eventlist = [1 2 3 4]
% 
% note: probabilities must 'sum to 1'

event = [];

n = length(eventlist);
m = length(probabilities);
if (n~=m), disp('Error! Size mismatch!'); return; end

intervals = zeros(n,1);

intervals(1) = probabilities(1); 
for i=2:n,
 intervals(i) = intervals(i-1) + probabilities(i);    
end


r     = rand;
i     = 1;
found = 0;
while (~found)
 if (r < intervals(i)), event = eventlist(i); found = 1; end;
 i = i + 1;
end

if (~found), disp('This violates axioms...!!!'); end;
end
