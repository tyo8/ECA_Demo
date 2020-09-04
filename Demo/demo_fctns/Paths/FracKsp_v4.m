% This function computes the ky- and kx- coordinates of a fractional Cart.
% sampling scheme. path.Ky and path.Kx are nRO x nPE double matrices with
% each column representing a single pulse(?) of the pulse sequence. e.g.
% Here, one column corresponds to one horizontal path across the k-space
% plane. The sampling density is uniform in the readout direction, and is
% designed to sample sparsely from high- and low-frequency spatial
% frequency bands in the phase-encode direction. Frequency bands should
% satisfy the Nyquist condition in their interiors. 
%
% If the gyro ratio and the acquisition parameters are provided, the
% function also computes the time stamps of the path coordinates.
% 
% The output is passed on to the k-space signal sampling module as an
% input.
% 
% For future versions, consider adding compatibility with pulse sequences.
%
%
% Version Control:
%   v1 adds a line (or two) at the center of k-space, but leaves total
%   acquisition time unchanged by removing equal numbers of lines from
%   other sampling segments, leaving gaps in the k-space weave
%
%   v2 adds several (2w+1) lines at the center of k-space without sacrificing
%   k-space coverage
%
%   v3 retains v1 k-space center behavior, but extends to 3 dimensions.
%
%   v4 simplifies Path structure to only contain the time-course data
%
% Inputs
%    AcqPars (struct)
%        TR (real scalar)
%        TE (real scalar)
%        FlipAngle (real scalar)
%        Resolution (real ndims-vector)
%        FOVsize (real ndims-vector)
%        nx/ny/nz (each is a real integer)
%        stddev (dims-sized array of noise variance)
%        onescantime
%        nscan
%        totalscantime
%    kernshift (real scalar) time from center to edge of kernel (ms)
%    ScanTime (real scalar) length of scan (in ms)
%    StartTime (real scalar, opt.) time after contrast injection (in ms)
%
% Outputs
%    path (struct)
%        Type (character) 'Cart'
%        Ky (double 2D array) ky-coordinates of Cartesian sampling path
%        Kx (double 2D array) kx-coordinates of Cartesian sampling path
%        Time (double 2D array) time-coordinates of Cartesian sampling path

function [Path,varargout] = FracKsp_v4(AcqPars, kernshift, varargin)
ScanTime = AcqPars.onescantime;

StartTime = 0;
if nargin > 2
    StartTime = varargin{1};
end

% Compute k-coordinates of Cartesian sampling path:
FOV = AcqPars.FOV;

% Compute K-coordinate grid of Cartesian sampling path:
kdims = FOV; %nPE_y, nRO, nPE_z

maxKx = 1/2;
Gyro = 42580; %gyromagnetic ratio (kHz/T)

TR = AcqPars.TR; %repetition time (milliseconds between each line acquisition)
TE = AcqPars.TE; %echo time (milliseconds from center to edge of k-space and back along readout dimension)
TS = maxKx / (Gyro * 10^(-5)); %time spent reading out each line (ms)
% TS = .5 / (Gyro * 10^(-5)); %time spent reading out each line (ms)
TDP = TE - TS / 2; %de-phasing time (ms)

tRO = TDP + TS * (1:kdims(2)) / kdims(2); %kdims(1) = nRO

% Compute time stamps according to aquisition parameters:
Path = zeros(kdims);
tPE_y = TR * (0:(kdims(1)-1)); %kdims(2) = nPE_y
tPE_z = TR * (0:(kdims(3)-1)); %kdims(3) = nPE_z; we assume that de-phasing in y-encoding also spoils z-phasing

for i=1:kdims(3)
    Path(:,:,i) = ones(kdims(1),1) * tRO + tPE_y' * ones(1,kdims(2)) + tPE_z(i) + TR*(kdims(2)-1)*(i-1);
end
path_length = Path(end)-StartTime;
scale_factor = AcqPars.onescantime/path_length;
Path = Path*scale_factor;

% Compute subdivision (and permutation) of Cartesian path that corresponds
% to the reconstruction kernels specified by KernPars

[nPE,nPart] = DefinePartitions(ScanTime,kernshift,FOV,TR);

base_idx = 1:nPE;
partition_num = 1 + mod(base_idx-1,nPart);
[~,fracidx] = sort(partition_num);

Path = Path + StartTime;

Path = permute(Path,[2 1 3]);
kdims_flipped = size(Path);
Path = reshape(Path(:,fracidx),kdims_flipped);
Path = permute(Path,[2 1 3]);

assert(isequal(kdims,FOV,size(Path)),...
    'mismatched dimensions in k-space sampling path')

[Path] = ReplicatePath(Path,AcqPars,StartTime);

if nargout > 2
    varargout{1} = fracidx;
end
end

function [nPE,nPart] = DefinePartitions(ScanTime,kernshift,FOV,TR)

FOV(2) = [];
real_scan_time = prod(FOV)*TR;

timescale_factor = real_scan_time/ScanTime;

ScanTime = ScanTime*timescale_factor;
kernshift = kernshift*timescale_factor;

nPE = floor(ScanTime/(TR)); %number of phase-encode lines acquired in entire scan
nPart = floor(ScanTime/(kernshift)); %number of partitions of a k-space sweep per reconstructed image
end

function [Path] = ReplicatePath(Path,AP,StartTime)

if isfield(AP,'nscan')
    repnum = AP.nscan;
elseif isfield(AP,'totalscantime') && isfield(AP,'onescantime')
    repnum = floor(AP.totalscantime/AP.onescantime);
else
    error('acquisition parameters do not specify total scan time')
end

p1 = Path;

for i=2:repnum
    p_new = p1 + max(Path(:)) - StartTime;
    Path(:,:,:,i) = p_new;
end
end