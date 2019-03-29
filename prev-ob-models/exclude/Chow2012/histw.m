function [n,bin] = histw(varargin)

[cax,args,nargs] = axescheck(varargin{:});

y = args{1};
x = args{2};
if nargs > 2
    w = args{3};
else
    w = ones(size(y));
end

bin = sort(x);
n = zeros(size(bin));
for i = 1:length(y)
    temp = find((bin-y(i))>=0);
    if length(temp) > 0
        n(temp(1)) = n(temp(1))+w(i);
    else
        N(end) = n(end)+w(i);
    end
end

end
