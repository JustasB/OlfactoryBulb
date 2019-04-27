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
