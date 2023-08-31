% function [ PSLR ] = PSLRindB( x0 )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
clc
a = get(gca,'Children');
xdata = get(a, 'XData');
ydata = get(a, 'YData');

for ii = 1:5
    beam = ydata{ii}; beam(end)
f = pc_freqwin( ydata{ii} );
pks = findpeaks(f);
pks = sort(pks);
PSLR = pks(length(pks))-pks(length(pks)-1);
PSLR_dB = pow2db(PSLR)
end
% end

