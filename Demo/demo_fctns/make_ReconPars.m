function [ReconPars] = make_ReconPars(PhantomPars,AcqPars)

ReconPars.convergethresh = 1e-8;
ReconPars.maxiter = 10^4;

ReconPars.smoothing = 1e-5;
ReconPars.weights = 1+PhantomPars.Function_Labels;

ReconPars.recontimeres = 250;
ReconPars.nimage = floor(AcqPars.totalscantime/ReconPars.recontimeres);
ReconPars.recontimeres = AcqPars.totalscantime/ReconPars.nimage;
end