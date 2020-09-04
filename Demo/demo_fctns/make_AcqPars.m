function [AcqPars] = make_AcqPars(Phant)
pdims = size(Phant.PhantomPars.Const_Image);

AcqPars.Resolution = ones(1,length(pdims));
AcqPars.nx = pdims(2);
AcqPars.ny = pdims(1);
AcqPars.nz = pdims(3);
AcqPars.FOV = pdims;

AcqPars.FlipAngle = 10*pi/180;

AcqPars.TR = 3.2;
AcqPars.TE = 1.6;

AcqPars.randomseed = rng;

AcqPars.acqtimeres = 50;

AcqPars.nscan = 12;
AcqPars.totalscantime = 42000;
AcqPars.onescantime = 3500;

AcqPars.phantomtimeres = 50;

AcqPars.stddev = Phant.PhantomPars.KspStdDev;
end