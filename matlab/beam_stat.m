function beam_stat(particles,scheduler)
% function beam_stat(particles)

%
% $Id: beam_stat.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
%
% Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
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

L=particles.rmax(particles.beamdir)-particles.rmin(particles.beamdir)-...
	particles.depth(particles.beamdir);
V=abs(particles.u(particles.beamdir))-...
	sqrt(2)*erfinv(particles.percent)*particles.vth(particles.beamdir);

tconst=L/V;

DIMR=[1:scheduler.dimr];
if tconst>0 & max(particles.t)>tconst,

	nmean=mean(particles.n(find(particles.t>tconst)));

	Ltotal=particles.rmax(particles.beamdir)-particles.rmin(particles.beamdir);
	V=abs(particles.u(particles.beamdir));
	corrseff=nmean/(abs(particles.flux(particles.beamdir))*Ltotal/V);

	ii=find(DIMR~=particles.beamdir);
	if scheduler.dimr==2,
		seff=particles.seff(ii);
	else
		seff=pi*prod(particles.seff(ii));
	end
	fprintf(1,'seff=%.2f\n', seff);
	fprintf(1,'seff corrected=%.2f\n', corrseff);
	ii=find(DIMR~=particles.beamdir);
	if scheduler.dimr==2,
		fprintf(1,'%% erf=%.3e\n', erf(corrseff/sqrt(2)/particles.width(ii(1))));
	else
		seff(1)=corrseff/pi/prod(particles.width(ii));
		seff(2)=corrseff/pi/prod(particles.width(ii));
		fprintf(1,'%% erf(%d)=%.3e\n',ii(1), ...
			erf(seff(1)/sqrt(2)/particles.width(ii(1))));
		fprintf(1,'%% erf(%d)=%.3e\n',ii(2), ...
			erf(seff(2)/sqrt(2)/particles.width(ii(2))));
	end

	plot(particles.t,particles.n,[tconst tconst], ...
		[0 nmean],[0 max(particles.t)],[nmean nmean]);

end
