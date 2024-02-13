function [frq,dfrq]=findfrq(t,y)
% function [frq,dfrq]=findfrq(t,y)

%
% $Id: findfrq.m,v 1.2 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2001-2011 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2.  of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

nt=length(t);

lt=max(t)-min(t);

nfft=1024;

frqs=[-nfft/2:nfft/2-1]'/(nfft-1)*(nt-1)/lt;

fty=fftshift(abs(fft(y-mean(y(:)),nfft)).^2);
subplot(111), plot(frqs,fty)

mxy=max(fty(:));
i=find(fty==mxy);

fprintf(1,'Calculated values for frq\n');
fprintf(1,'frq=%+.6f \n', frqs(i));
fprintf(1,'dfrq=%f\n', max(diff(frqs)));

frq=frqs(i);

dfrq=max(diff(frqs));
