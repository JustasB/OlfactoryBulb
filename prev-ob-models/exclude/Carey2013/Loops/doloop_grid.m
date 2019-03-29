function doloop_grid(tracefilepath, modelfnc, idx, varargin)
% DOLOOP_GRID runs the selected juxtaglomerular model for each ORN input trace
% (glomerulus) in the selected input file, modifying model parameters as
% specified across a 2D grid, then saves the model output to a file.
%
%   INPUTS: 
%           tracefilepath: fullpath to ORN input data
%           modelfnc: name of juxtaglomerular model file
%           idx: unused except in naming output file
%           varargin: list of parameters/values to change from defaults
%           when running the model. Should be given as (par1, grid1, par2,
%           grid2), where gridN is the vector of values to be used for parN

    [~,tracename,~] = fileparts(tracefilepath);
    modelname = func2str(modelfnc);
    
    paroverrides = varargin;
    params = paroverrides(1:2:length(paroverrides));
    grids = paroverrides(2:2:length(paroverrides));
    nConds = cellfun(@length,grids);

    S = load(tracefilepath);
    
    %special check for "trialtraces" variable - for awake orn input data
    if isfield(S,'trialtraces')
        [isTrialSelected, i] = ismember('trialindex', params);
        if isTrialSelected
            S.traces = S.trialtraces{grids{i}};
            trialname = S.trialnames{grids{i}};
            params(i) = [];
            grids(i) = [];
        end
    else
        isTrialSelected = false;
    end
    
    nGlom = size(S.traces,1);
    
    condition_index = cell(1,length(grids)+1);
    [condition_index{:}] = ind2sub([nGlom nConds],1:(prod(nConds)*nGlom));
    conditions = cell2mat(condition_index');
    glomIndices = conditions(1,:);
    
    paroverridegrid = cell(size(conditions,2),length(params)*2);
    for n = 1:length(params)
        paroverridegrid(:,2*n-1) = params(n);
        paroverridegrid(:,2*n) = num2cell(grids{n}(conditions(n+1,:)));
    end
    
    total_start = tic;
    
    for i = 1:size(paroverridegrid,1)
        disp(['Integrating condition ' num2str(i) '...'])
        data(i) = modelfnc(S.traces(glomIndices(i),:), S.samplingrate, paroverridegrid{i,:});
        disp(['... finished condition ' num2str(i) '.'])
    end
    
    %data = rmfield(data, {'T','X'});
    
    c = num2cell(glomIndices);
    [data.glomTrace] = deal(c{:});
    extrastr = '';
    for n = 1:2:size(paroverridegrid,2)
        [data.(paroverridegrid{1,n})] = deal(paroverridegrid{:,n+1});
        if length(paroverrides{n+1}) == 1
            if iscell(paroverrides{n+1})
                st = paroverrides{n+1}{1};
            else
                st = num2str(paroverrides{n+1});
            end
            extrastr = [extrastr '__' paroverrides{n} regexprep(st,'\.','_')];
        end
    end
    
    if isTrialSelected
        extrastr = [trialname '__' extrastr];
    end
    
    disp('Total time ...');
    toc(total_start)
    
    save(sprintf('%s_%s%s_grid_%u',modelname,tracename,extrastr,idx),'data')
end