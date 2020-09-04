% t is parameter time and should be given in units of milliseconds; however, time
% units in the AIF function parameters are evaluated in minutes, so 't' is
% re-scaled within the function (t = t/60).

function [con,maxcon] = AIF(t,varargin)
t = (t/6e4).*(t>0); %convert from milliseconds to minutes

popG1p = {.809,.0563,.1705};
popG2p = {.33, .132,.365};
popSp  = {1.05, .1685, 38.078, 0.483};
def_norm = false;

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'Sp',popSp)
addParameter(p,'G1p',popG1p)
addParameter(p,'G2p',popG2p)
addParameter(p,'norm',def_norm)

parse(p,varargin{:})

Sp = p.Results.Sp;
G1p = p.Results.G1p;
G2p = p.Results.G2p;
norm = p.Results.norm;

con = Sp{1}.*exp(-Sp{2}.*t)./(1+exp(-Sp{3}.*(t-Sp{4})))...
    + G1p{1}./(G1p{2}.*sqrt(2*pi)) * exp(-(t-G1p{3}).^2./(2*G1p{2}.^2))...
    + G2p{1}./(G2p{2}.*sqrt(2*pi)) * exp(-(t-G2p{3}).^2./(2*G2p{2}.^2));
        
maxcon = Sp{1}*exp(-Sp{2}*G1p{3})./(1+exp(-Sp{3}*(G1p{3}-Sp{4})))...
    + G1p{1}/(G1p{2}*sqrt(2*pi)) * exp(-(G1p{3}-G1p{3}).^2/(2*G1p{2}^2))...
    + G2p{1}/(G2p{2}*sqrt(2*pi)) * exp(-(G1p{3}-G2p{3}).^2/(2*G2p{2}^2));

if norm
    con = con/maxcon;
end

end