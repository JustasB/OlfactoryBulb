function val = focality(P, metric)
    [Nc Ns] = size(P);
    tempval = [];
    for sti = 1:Ns
        MCact = rec(P(:,sti),0);
        MCacttemp = MCact;
        MCacttemp(find(MCact<=max(MCact)/2)) = 0;
        nor_down = 0;
        down = 0;
        for ii = 1:length(MCact)
            for jj = ii:length(MCact)
                down = down + metric(ii,jj);
                nor_down = nor_down+1;
            end
        end
        down = down/nor_down;
        nor_up = 0;
        up = 0;
        for ii = 1:length(MCact)
            for jj = ii:length(MCact)
                up = up + MCacttemp(ii)*MCacttemp(jj)*metric(ii,jj);
                nor_up = nor_up + MCacttemp(ii)*MCacttemp(jj);
            end
        end
        up = up/nor_up;
        tempval = [tempval 1-up/down];
    end
    val = tempval;
end % focality