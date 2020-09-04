% Runs a demo of the virtual scanner and time-tagged reconstruction process
% on a toy patient phantom.
% 
% Inputs:
%     - pathtype, an integer in the set {0,1,2}; 0 is a standard Cartesian
%       acquisition, 1 is an inter-aliased Cartesian acquisition, and 2 is
%       a random-order Cartesian acquisition (all paths are 3D sequences)
% 
% Outputs:
%     - Img_rec, the reconstructed image
%     - resnorms, the gradient norm at each iteration
% 
% Displays:
%     - A comparison sample slice from gold standard, TTR-reconstructed,
%       and IFFT-reconstructed images, at each full scan point
%     - A plot of splined signal time-course data from a few sample voxels
%       from gold-standard, TTR-reconstructed, and IFFT-reconstructed
%       images

function [Img_rec] = run_demo(pathtype)
start = datetime;
addpath(genpath('demo_fctns'))

switch pathtype
    case 0
        pathname = 'Standard Cartesian';
    case 1
        pathname = 'Inter-aliased Cartesian';
    case 2
        pathname = 'Random Cartesian';
end

%in demo, runs from install site
dirname = 'Demo';
phantfname = 'demo_phant.mat';
phantname = 'Phant';

pfname = fullfile(dirname,phantfname);
vars = load(pfname,phantname);
Phant = vars.(phantname);
clear vars

disp(newline)
disp(['Welcome to the ' pathname ' demo!'])

disp(newline)
disp('Generating acquisition parameters...')
AcqPars = make_AcqPars(Phant);
disp('Done.')

disp(newline)
disp('Generating K-space path...')
Path = make_Path(AcqPars,pathtype);
disp('Done.')

disp(newline)
disp('Generating K-space signal...')
tic
Y = GenerateKspaceSignal(Phant.PhantomEvalFn,Phant.PhantomPars,Path,AcqPars);
disp(['Done in ' num2str(toc) ' seconds.'])

ReconPars = make_ReconPars(Phant.PhantomPars,AcqPars);

disp(newline)
disp('Computing Time-Tagged Reconstruction...')
tic
Img_rec = ReconAlg(Y,Path,AcqPars,ReconPars);
disp(['Done in ' num2str(toc) ' seconds.'])
disp(newline)


%%compare and display test results
do_and_show_analysis(Phant,Y,Img_rec,Path,AcqPars,ReconPars,pathtype);

finish = datetime;
elapsed = etime(datevec(finish),datevec(start));

disp(['This demo ran in only ' num2str(elapsed) ' seconds!'])
end

function [Path] = make_Path(AcqPars,pathtype)
switch pathtype
    case 0
        Path = PathCart_v5(AcqPars);
    case 1
        Path = FracKsp_v4(AcqPars,250);
    case 2
        Path = RandCart_v2(AcqPars);
end
end

function [] = do_and_show_analysis(Phant,Y,Img_ttr,Path,AcqPars,ReconPars,pathtype)
Img_ttr = real(Img_ttr);

plot_enh_curves(Phant,Y,Img_ttr,Path,AcqPars,ReconPars)

display_imgs(Phant,Y,Img_ttr,Path,AcqPars,ReconPars,pathtype)

end

function [] = plot_enh_curves(Phant,Y,Img_ttr,Path,AcqPars,ReconPars)
vox_const = [10, 13, 10];
vox_AIF = [8, 2, 7];
vox_EMM = [9, 14, 9];

dts = AcqPars.onescantime;
t_ifft = (min(Path(:)):dts:(dts*size(Y,4))) + dts/2;
Img_ifft = real(ifft3d(Y));

dtr = ReconPars.recontimeres;
t_ttr = (min(Path(:)):dtr:(dtr*size(Img_ttr,4))) + dtr/2;

t_gsi = t_ttr;
TrueImg = zeros(size(Img_ttr));
for i = 1:length(t_ttr)
    TrueImg(:,:,:,i) = Phant.PhantomEvalFn(Phant.PhantomPars,AcqPars,t_ttr(i));
end

const_ttr_i = squeeze(Img_ttr(vox_const(1),vox_const(2),vox_const(3),:));
const_ttr = spline(t_ttr,const_ttr_i,t_gsi);
const_ifft_i = squeeze(Img_ifft(vox_const(1),vox_const(2),vox_const(3),:));
const_ifft = spline(t_ifft,const_ifft_i,t_gsi);
const_gsi = squeeze(TrueImg(vox_const(1),vox_const(2),vox_const(3),:));

figure,plot(t_gsi,const_gsi)
hold on
plot(t_ttr,const_ttr_i,'-r.')
plot(t_ifft,const_ifft_i,'-.mo');
hold on
title('Voxel with constant value')
xlabel('Time (ms)')
legend('Ground Truth','TT Recon','IFFT Recon')

AIF_ttr_i = squeeze(Img_ttr(vox_AIF(1),vox_AIF(2),vox_AIF(3),:));
AIF_ttr = spline(t_ttr,AIF_ttr_i,t_gsi);
AIF_ifft_i = squeeze(Img_ifft(vox_AIF(1),vox_AIF(2),vox_AIF(3),:));
AIF_ifft = spline(t_ifft,AIF_ifft_i,t_gsi);
AIF_gsi = squeeze(TrueImg(vox_AIF(1),vox_AIF(2),vox_AIF(3),:));

figure,plot(t_gsi,AIF_gsi)
hold on
plot(t_ttr,AIF_ttr_i,'-r.')
plot(t_ifft,AIF_ifft_i,'-.mo');
hold off
title('AIF Voxel')
xlabel('Time (ms)')
legend('Ground Truth','TT Recon','IFFT Recon')

EMM_ttr_i = squeeze(Img_ttr(vox_EMM(1),vox_EMM(2),vox_EMM(3),:));
EMM_ttr = spline(t_ttr,EMM_ttr_i,t_gsi);
EMM_ifft_i = squeeze(Img_ifft(vox_EMM(1),vox_EMM(2),vox_EMM(3),:));
EMM_ifft = spline(t_ifft,EMM_ifft_i,t_gsi);
EMM_gsi = squeeze(TrueImg(vox_EMM(1),vox_EMM(2),vox_EMM(3),:));

figure,plot(t_gsi,EMM_gsi)
hold on
plot(t_ttr,EMM_ttr_i,'-r.')
plot(t_ifft,EMM_ifft_i,'-.mo');
hold off
title('EMM Voxel')
xlabel('Time (ms)')
legend('Ground Truth','TT Recon','IFFT Recon')
end

function [] = display_imgs(Phant,Y,Img_ttr_all,Path,AcqPars,ReconPars,pathtype)
dts = AcqPars.onescantime;
t_ifft = (min(Path(:)):dts:(dts*size(Y,4))) + dts/2;
Img_ifft = real(ifft3d(Y));

dtr = ReconPars.recontimeres;
t_ttr_all = (min(Path(:)):dtr:(dtr*size(Img_ttr_all,4))) + dtr/2;
for i = 1:length(t_ifft)
    [~,t_ttr_idx(i)] = min(abs(t_ttr_all-t_ifft(i)));
end
Img_ttr = Img_ttr_all(:,:,:,t_ttr_idx);

TrueImg = zeros(size(Y));
for i = 1:length(t_ifft)
    TrueImg(:,:,:,i) = Phant.PhantomEvalFn(Phant.PhantomPars,AcqPars,t_ifft(i));
end

zslice = round(size(Img_ifft,3)/2);
nt = size(Img_ifft,4);

True = squeeze(TrueImg(:,:,zslice,:));
IFFT = squeeze(Img_ifft(:,:,zslice,:));
TTR = squeeze(Img_ttr(:,:,zslice,:));

Img_composite = cat(3,True,IFFT,TTR);
maxsig = max(abs(Img_composite(:)));

switch pathtype
    case 0
        pathname = 'standard Cartesian';
    case 1
        pathname = 'inter-aliased Cartesian';
    case 2
        pathname = 'random Cartesian';
end

figure, montage(Img_composite/maxsig,'size',[3 nt])
title(['Reconstructed from ' pathname ' sampling'])
end