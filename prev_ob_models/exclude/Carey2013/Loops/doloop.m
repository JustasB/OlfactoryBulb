function doloop(tracefilepath, ORNscale, modelfnc, varargin)
% DOLOOP runs the selected juxtaglomerular model for each ORN input trace
% (glomerulus) in the selected input file, modifying model parameters as
% specified, then saves the model output to a file.
%
%   INPUTS: 
%           tracefilepath: fullpath to ORN input data
%           ORNscale: ORN gain premultiplying the ORN input data
%           modelfnc: name of juxtaglomerular model file
%           varargin: list of parameters/values to change from defaults
%           when running the model. Should be given as par1, val1, par2,
%           val2, ...

    [~,tracename,~] = fileparts(tracefilepath);
    modelname = func2str(modelfnc);
    scalingstr = regexprep(num2str(ORNscale),'\.','_');
    
    S = load(tracefilepath);
    
    nGlom = size(S.traces,1);
    
    total_start = tic;
    
    for i = 1:nGlom
        disp(['Integrating glomerulus ' num2str(i) '...'])
        data(i) = modelfnc(ORNscale * S.traces(i,:), S.samplingrate, varargin{:});
        disp(['... finished glomerulus ' num2str(i) '.'])
    end
    
    %data = rmfield(data, {'T','X'});
    [data.ORNscale] = deal(ORNscale);
    paroverrides = varargin;
    extrastr = '';
    for p = 1:2:length(paroverrides)
        [data.(paroverrides{p})] = deal(paroverrides{p+1});
        extrastr = [extrastr '_' paroverrides{p} regexprep(num2str(paroverrides{p+1}),'\.','_')];
    end
    
    disp('Total time ...');
    toc(total_start)
    
    save(sprintf('%s_%s_gain%s%s',modelname,tracename,scalingstr,extrastr),'data')
end
