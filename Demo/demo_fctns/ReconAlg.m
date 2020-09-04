function [EstSignal,resnorms] = ReconAlg(KspaceSignal,Path,AcqPars,ReconPars)

%% 
% Setup
% renormalize data (smaller voxel values converge slightly faster)
maxsig = max(abs(KspaceSignal(:)));
KspaceSignal = KspaceSignal/maxsig;

% smoothing matrix

Dmat = spdiags([ReconPars.smoothing + ...
    [1;5;6*ones(ReconPars.nimage-4,1);5;1],...
    [0;-2;-4*ones(ReconPars.nimage-3,1);-2],...
    [-2;-4*ones(ReconPars.nimage-3,1);-2;0],...
    ones(ReconPars.nimage,2)],[0,1,-1,2,-2],...
    ReconPars.nimage,ReconPars.nimage);

% mask of observed & unobserved times

ObsMask = false(AcqPars.nx,AcqPars.ny,AcqPars.nz,ReconPars.nimage);
for iimage = 1:ReconPars.nimage
    timerange = (AcqPars.totalscantime/ReconPars.nimage)*[iimage-1,iimage];
    if(iimage==1)
        timerange(1) = -inf;
    end
    if(iimage==ReconPars.nimage)
        timerange(2) = inf;
    end
    ObsMask(:,:,:,iimage) = any((Path > timerange(1)) & ...
        (Path <= timerange(2)),4);
end

% rearrange KspaceSignal to the new time resolution
KspaceSignal_obs = zeros(size(ObsMask));
for ix = 1:AcqPars.nx
    for iy = 1:AcqPars.ny
        for iz = 1:AcqPars.nz
            KspaceSignal_obs(ix,iy,iz,ObsMask(ix,iy,iz,:)) = ...
                KspaceSignal(ix,iy,iz,:);
        end
    end
end

% weights stretched over time axis
weights_stretch = repmat(ReconPars.weights,1,1,1,ReconPars.nimage);

%%
% Reconstructing the signal

% Solve (Dinv x (F*Winv*F'))_{obs,obs} * v = KspaceSignal
% Then set EstSignal = (Dinv x (Winv*F'))_{...,obs} * v

% Conjugate gradient descent to solve the linear system for v
mult_by_matrix = @(v) ...
    reshape(ObsMask.*(fft3d((ifft3d(...
    reshape(reshape(ObsMask(:).*v,[],ReconPars.nimage)/Dmat,...
    AcqPars.nx,AcqPars.ny,AcqPars.nz,[])))./weights_stretch)),[],1);

[v,flag,relres,iter,resnorms] = pcg(mult_by_matrix,KspaceSignal_obs(:),...
    ReconPars.convergethresh,ReconPars.maxiter);

% Compute EstSignal as a function of v
v_times_Dinv = reshape(v,[],ReconPars.nimage)/Dmat;
EstSignal = ifft3d(reshape(v_times_Dinv,...
    AcqPars.nx,AcqPars.ny,AcqPars.nz,[]))./weights_stretch;
    
% rescale voxel values
EstSignal = EstSignal*maxsig;
end