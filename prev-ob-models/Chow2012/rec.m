function out = rec(V,r,smooth)
    if nargin <3
        smooth = 0.1;
    end
    
    if smooth == 0
        out = (abs(V-r)+V-r)/2;
    else
        x = V-r;
        y = zeros(size(x));
        b = 1/smooth;
        c = log(2*b)/b;

        IX = find(x<=0);
        y(IX) = exp(b*(x(IX)-c));
        IX = find(x>0);
        y(IX) = x(IX) + exp(-b*(x(IX)+c));
        out = y;
    end
end % rec