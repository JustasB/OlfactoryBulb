function S = get_odor(name)

    odor_raw = imread(strcat('odor/',name,'.png'));
    odor = zeros(size(odor_raw(:,:,1)));

    for i = 1:size(odor,1)
        for j = 1:size(odor,2)
            if odor_raw(i,j,1) == 255
                if odor_raw(i,j,2) == 255
                    if odor_raw(i,j,3) == 255
                        odor(i,j) = -2;
                    end
                end
            end
            if odor_raw(i,j,1) == 0
                if odor_raw(i,j,2) == 0
                    if odor_raw(i,j,3) == 255
                        odor(i,j) = -1;
                    end
                end
            end
            if odor_raw(i,j,1) == 0
                if odor_raw(i,j,2) == 255
                    if odor_raw(i,j,3) == 255
                        odor(i,j) = 0;
                    end
                end
            end
            if odor_raw(i,j,1) == 0
                if odor_raw(i,j,2) == 193
                    if odor_raw(i,j,3) == 0
                        odor(i,j) = 1;
                    end
                end
            end
            if odor_raw(i,j,1) == 255
                if odor_raw(i,j,2) == 255
                    if odor_raw(i,j,3) == 0
                        odor(i,j) = 2;
                    end
                end
            end
            if odor_raw(i,j,1) == 255
                if odor_raw(i,j,2) == 133
                    if odor_raw(i,j,3) == 0
                        odor(i,j) = 3;
                    end
                end
            end
            if odor_raw(i,j,1) == 255
                if odor_raw(i,j,2) == 0
                    if odor_raw(i,j,3) == 0
                        odor(i,j) = 4;
                    end
                end
            end
            if odor_raw(i,j,1) == 0
                if odor_raw(i,j,2) == 0
                    if odor_raw(i,j,3) == 0
                        odor(i,j) = 5;
                    end
                end
            end
        end
    end
    
    odor(1,:) = NaN;
    odor(:,1) = NaN;
    for i = 2:size(odor,1)
        for j = 2:size(odor,2)
            if odor(i,j) == -2
                if isnan(odor(i-1,j))
                    odor(i,j) = NaN;
                elseif isnan(odor(i,j-1))
                    odor(i,j) = NaN;
                end
            end
        end
    end
    odor(end,:) = NaN;
    odor(:,end) = NaN;
    for i = size(odor,1)-1:-1:1
        for j = size(odor,2)-1:-1:1
            if odor(i,j) == -2
                if isnan(odor(i+1,j))
                    odor(i,j) = NaN;
                elseif isnan(odor(i,j+1))
                    odor(i,j) = NaN;
                end
            end
        end
    end
    
    S = odor;

end