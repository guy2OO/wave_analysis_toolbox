function [waveLength] = Wave_Length(period,depth)

%% Inputs
%
%   period - wave period (s)
%
%   depth - water depth (m)
%

%% Initialize known variables and file data
g = 9.807;

% Calculate angular frequency omega
omega = 2 * pi./period;

% Calculate approximate wave length
ApxWaveLengthD = (g/(2*pi)).*period.^2;
ApxWaveLengthS = sqrt(g.*depth.*period);

% Solve for wave length using omega and gravity
% syms L
% waveLength(i) = vpasolve(omega^2 == g*(2*pi/L)*tanh((2*pi*depth)/L),L);
% 
% waveLength(i) = abs(waveLength(i));

% Create annonymos function equal to zero to solve for with fsolve and iterate
waveLength = zeros(1,length(period));
for i = 1:length(period)
    dipf=@(L) omega(i).^2 -g.*(2*pi./L).*tanh((2*pi.*depth)./L);

    waveLength(i) = fsolve(dipf,min([ApxWaveLengthD(i) ApxWaveLengthS(i)]));
end