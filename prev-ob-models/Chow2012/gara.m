function [S, coord, metric, name] = gara(Ns,MCperGlom,name_list,sim,range)
    
    % read raw pattern
    name_list = strrep(name_list,' ','#');
    rem = name_list;
    [token, rem] = strtok(rem,'#');
    temp = get_odor(token);
    odor = zeros(Ns,size(temp,1),size(temp,2));
    odor(1,:,:) = temp;
    i = 1;
    while length(rem)>1
        i = i+1;
        [token, rem] = strtok(rem,'#');
        temp = get_odor(token);
        odor(i,:,:) = temp;
    end
    raw_Ns = size(odor,1);
    name_list = strrep(name_list,'(_)','(-)');
    name_list = strrep(name_list,'_',' ');
    rem = name_list;
    [token,rem] = strtok(rem,'#');
    max_length = length(token);
    for i = 2:raw_Ns
        [token,rem] = strtok(rem,'#');
        if length(token)>max_length
            max_length = length(token);
        end
    end
    rem = name_list;
    name = cell(1,raw_Ns);
    for i = 1:raw_Ns
        [token, rem] = strtok(rem,'#');
        temp = token;
        name{i} = temp;
    end
    
%     % plot original image
%     figure(1000); clf;
%     for i = 1:Ns
%         subplot(3,4,i);
%         imagesc(reshape(odor(i,:,:),size(odor,2),size(odor,3)));
%         caxis([min(min(min(odor))) max(max(max(odor)))]);
%         axis off; axis image;
%         title(name{i});
%     end

    % chop image
    if nargin > 4
        odor = odor(:,range(1):range(2), range(3):range(4));
%         for i = 1:raw_Ns
%             subplot(2,8,8+i);
%             imagesc(reshape(odor(i,:,:),size(odor,2),size(odor,3)));
%             caxis([min(min(min(odor))) max(max(max(odor)))]);
%             axis off; axis image;
%         end
    end
    
    % reduce pix
    for sti = 1:raw_Ns
        odorsti = reshape(odor(sti,:,:),size(odor,2),size(odor,3));
        for i = 1:sim
            tempi = odorsti(i:sim:sim*floor(end/sim), i:sim:sim*floor(end/sim));
            if i==1
                temp = tempi;
            else
                temp = max(temp, tempi);
            end
        end
        if sti == 1
            odor_redu = zeros(raw_Ns,size(temp,1),size(temp,2));
            odor_redu(1,:,:) = temp;
        else
            odor_redu(sti,:,:) = temp;
        end
    end
    
    % find revalent pix
    odor_sum = reshape(sum(odor_redu,1),size(odor_redu,2),size(odor_redu,3));
    odor_sum_1d = reshape(odor_sum, 1, size(odor_sum,1)*size(odor_sum,2));
    IX = find(odor_sum_1d>-3);
    odor_sum_1d = odor_sum_1d(IX);
    
    % cal coord1mc
    coord1mc1 = (ones(size(odor_sum,2),1)*(size(odor_sum,1):-1:1))';
    coord1mc2 = ((1:size(odor_sum,2))'*ones(1,size(odor_sum,1)))';
    coord1mc1 = reshape(coord1mc1,1,size(coord1mc1,1)*size(coord1mc1,2));
    coord1mc2 = reshape(coord1mc2,1,size(coord1mc2,1)*size(coord1mc2,2));
    coord1mc = [coord1mc2; coord1mc1];
    coord1mc = coord1mc(:,IX);
    
    % select relavent pix for each odor
    for sti = 1:raw_Ns
        if sti == 1
            odor_redu1d = zeros(raw_Ns,length(IX));
        end
        temp = reshape(odor_redu(sti,:,:),1,size(odor_sum,1)*size(odor_sum,2));
        odor_redu1d(sti,:) = temp(IX);
    end
    
    % sort according to strength
    % [temp,IX] = sort(sum(odor_redu1d));
    [temp,IX] = sort(mean(odor_redu1d(1:2,:))-mean(odor_redu1d(3:4,:))+0.01*sum(odor_redu1d));
    odor_redu1d = odor_redu1d(:,IX(end:-1:1));
    coord1mc = coord1mc(:,IX(end:-1:1));
    
    % metric
    metric1mc = zeros(size(odor_redu1d,2),size(odor_redu1d,2));
    for i = 1:size(odor_redu1d,2)
        for j = 1:size(odor_redu1d,2)
            metric1mc(i,j) = norm(coord1mc(:,i)-coord1mc(:,j));
        end
    end
 
    % return
    S1mc = odor_redu1d'/max(max(odor_redu1d));
    S1mc = S1mc(:,1:Ns);
    
    % mult. MC per glom
    coord = zeros(2,MCperGlom*size(coord1mc,2));
    metric = zeros(MCperGlom*size(metric1mc,1));
    S = zeros(MCperGlom*size(S1mc,1),size(S1mc,2));
    for i = 1:MCperGlom
        coord(:,i:MCperGlom:end) = coord1mc;
        S(i:MCperGlom:end,:) = S1mc;
        metric(i:MCperGlom:end, 1:MCperGlom:end) = metric1mc;
    end
    for i = 2:MCperGlom
        metric(:, i:MCperGlom:end) = metric(:, 1:MCperGlom:end);
    end
end