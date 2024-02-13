function profile=hydrostatic(ng,rmin,rmax,m,G,T)
% function profile=hydrostatic(ng,rmin,rmax,m,G,T)
% 
% For instance to generate the profile to be used in share/gravity_pulse2d.conf
%
%   profile=hydrostatic([145,337],[0,0],[95,250],1,[0,-1e-2],10);
%   write_nhbgrd_model(profile,'matlab/hydro');

%
% $Id: hydrostatic.m,v 1.6 2018/03/08 15:27:56 patrick Exp $
%
% Copyright (c) 2016 Patrick Guio <patrick.guio@gmail.com>
% All Rights Reserved.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.


profile = zeros(ng);


dir = find(G~=0);
if isempty(dir), % If G=0 assume direction to be y
  dir = 2;
end

H = m*G(dir)/T;


a = rmin(dir);
b = rmax(dir);

rs = linspace(rmin(dir), rmax(dir), ng(dir));

% normalised so that 1/(b-a)\int_a^b ns dr = 1;
ns = H*(b-a)/(exp(H*b)-exp(H*a))*exp(H*rs);
if ~isfinite(ns),
  ns = ones(size(rs));
end

fprintf(1,'1/(b-a)\\int_a^b ns dr = %.6f\n', ...
  sum((ns(1:end-1)+ns(2:end)).*diff(rs)/2 )/(b-a));

if 1,
  dr = rmax-rmin;
  rs0 = linspace(0, dr(dir), ng(dir));
  nrmin = -dr(dir)*H/(1-exp(dr(dir)*H));
	if ~isfinite(nrmin), % dealing with G=0
	  nrmin = 1; 
	end;
  nrmax = nrmin*exp(dr(dir)*H);
	if ~isfinite(nrmax), % dealing with G=0
	  nrmax = 1; 
	end;
  ns0 = nrmin*exp(H*rs0);
	if G(dir)>0, 
	  fprintf(1,'nmin=%.6f, nmax=%.6f lmin=%d\n',...
	          nrmin,nrmax,fix(1/nrmin)); 
	else, 
	  fprintf(1,'nmin=%.6f, nmax=%.6f lmin=%d\n',...
		        nrmax,nrmin,fix(1/nrmax)); 
	end
  fprintf(1,'1/(b-a)\\int_a^b ns dr = %.6f\n', ...
    sum((ns0(1:end-1)+ns0(2:end)).*diff(rs0)/2 )/(b-a));
  plot(ns,rs,ns0,rs0+rmin(dir))
	xlabel('n(z)'); ylabel('z');
	fprintf(1,'max|ns-ns0| = %.6e\n', max(abs(ns-ns0)))
end

switch(length(ng))

  case 2,
	   switch dir
		   case 1, profile = repmat(reshape(ns,[ng(1),1]),[1,ng(2)]);
		   case 2, profile = repmat(reshape(ns,[1,ng(2)]),[ng(1),1]);
     end
	case 3,
	   switch dir
		   case 1, profile = repmat(reshape(ns,[ng(1),1,1]), [1,ng(2:3)]);
		   case 2, profile = repmat(reshape(ns,[1,ng(2),1]), [ng(1),1,ng(3)]);
		   case 3, profile = repmat(reshape(ns,[1,1,ng(3)]), [ng(1:2),1]);
     end

end
 
