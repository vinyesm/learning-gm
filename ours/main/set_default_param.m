function param = set_default_param(param)
%% set default values to param 

if ~isfield(param,'k')
    param.k = 0;
end

if ~isfield(param,'q')
    param.q = param.k;
end


if ~isfield(param,'PSD')
    param.PSD = true;
end

if ~isfield(param,'nbMainLoop')
    param.nbMainLoop = 100;
end


if ~isfield(param,'powerIter')
    param.powerIter = 100;
end
if ~isfield(param,'stPtPowerIter')
    param.stPtPowerIter = 100;
end

if ~isfield(param,'epsStop')
    param.epsStop = 1e-6;
end

if ~isfield(param,'innerLoopIter')
    param.innerLoopIter = 100;
end

if ~isfield(param,'niterPS')
    param.niterPS =  5000;
end

if ~isfield(param,'PSdualityEpsilon')
    param.PSdualityEpsilon = 1e-6;
end

if ~isfield(param,'verbose')
    param.verbose = 0;
end

if ~isfield(param,'debug')
    param.debug = 0;
end

if ~isfield(param,'sloppy')
    param.sloppy = 0;
end

if ~isfield(param,'Sfixed')
    param.Sfixed = false;
end

if ~isfield(param,'Sstar')
    param.Sstar = 0;
end

if ~isfield(param,'no_l1')
    param.no_l1 = false;
end

end