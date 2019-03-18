function val = corr(M,type,S)
    if nargin < 2
        type = 1;
    end
    if type == 1
        val = corrcoef(M);
    elseif type == 2
        out = eye(size(M,2));
        M2 = zeros(size(M));
        for i = 1:size(M,2)
            if norm(M(:,i)) == 0
                out = NaN*out;
                return;
            end
            M2(:,i) = M(:,i)/norm(M(:,i));
        end
        val = M2'*M2;
    end
end