function [con,maxcon] = EMM(t,varargin)
t = t/1000; %convert from milliseconds to seconds
sampEMMp = {1.066471267, 0.168860905, 6.10530657};

p = inputParser;
p.KeepUnmatched = true;
addParameter(p,'EMMp',sampEMMp)

parse(p,varargin{:})
EMMp = p.Results.EMMp;

con = EMMp{1}.*(1-exp(-EMMp{2}.*(t-EMMp{3}))).*(t > EMMp{3});
maxcon = max(EMMp{1}(:));
end