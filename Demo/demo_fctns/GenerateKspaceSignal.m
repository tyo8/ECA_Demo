function [KspaceSignal, AcqPars] = ...
    GenerateKspaceSignal(PhantomEvalFn,PhantomPars,Path,AcqPars)

KspaceSignal = zeros(AcqPars.nx,AcqPars.ny,AcqPars.nz,AcqPars.nscan);

% perform Fourier transform at each measurement time,
% and extract the measurement at the corresponding point in k-space


ntime = round(AcqPars.totalscantime / AcqPars.phantomtimeres);
timegap = AcqPars.totalscantime / ntime;

kspace_pre = fft3d(PhantomEvalFn(PhantomPars,AcqPars,0));
kspace_post = fft3d(PhantomEvalFn(PhantomPars,AcqPars,timegap));

for itime = 1:ntime
    path_mask = logical((Path >= timegap*(itime-1)) & (Path <= timegap*itime));
    linear_combo_coef = (timegap*itime - Path) / timegap;
    KspaceSignal = KspaceSignal + path_mask .* (kspace_pre .* linear_combo_coef + kspace_post .* (1-linear_combo_coef));
    if(itime < ntime)
        kspace_pre = kspace_post;
        kspace_post = fft3d(PhantomEvalFn(PhantomPars,AcqPars,timegap*(itime+1)));
    end
end


% add complex Gaussian noise
if isfield(AcqPars,'randomseed')
    rng(AcqPars.randomseed)
else
    AcqPars.randomseed = rng;
end
KspaceSignal = KspaceSignal + ...
    randn(size(KspaceSignal)).*repmat(AcqPars.stddev,1,1,1,AcqPars.nscan);
end