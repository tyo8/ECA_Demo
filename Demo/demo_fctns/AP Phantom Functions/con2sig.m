% This function converts a generic contrast concentration measurement to a
% signal measurement. Input 'con' (contrast concentration) is given in
% units of mmol

function [sig] = con2sig(con,varargin)

p = inputParser;

addParameter(p,'FA',10*pi/180)
addParameter(p,'TR',3.2)
addParameter(p,'TE',1.6)

parse(p,varargin{:})

FA = p.Results.FA; %flip angle (rad)
TR = p.Results.TR; %repetition time (ms)

if FA > pi
    FA = FA*pi/180; %ensures FA is expressed in radians and not degrees
end

r1 = 3.4/1000; %(mmol*ms)^-1
R10 = 1/1.8/1000; %inverse native T1 (ms^-1)
R1 = R10 + r1*con; %ms^-1
T1 = 1./R1; %T1 relaxation time (ms)
% T1(R1 == 0) = 0;

S_c = exp(-TR./T1);
% S_c(T1 == 0) = 0;

T2_st = 250; %relaxation time (ms)
TE = 1.6; %echo time (ms)
p = 1;

sig = p*exp(-TE/T2_st)*sin(FA)*(1 - S_c)./(1 - cos(FA)*S_c).*(con > 0);
end