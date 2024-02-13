function energy3d
% function energy3d

%
% $Id: energy3d.m,v 1.2 2011/03/26 15:36:04 patrick Exp $
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

files=str2mat('pb0m0_3d.hdf','pb0m1_3d.hdf','ib0_3d.hdf','fb0_3d.hdf');

kei=[];
ese=[];
for i=1:size(files,1),

	filename = deblank(files(i,:));
	if exist(filename,'file')==2,
		filename=which(filename);
	else
		error(sprintf('Cannot find file %s', filename));
	end

	[tt,t,fx]=readvector(filename, 'ese');
	ese=[ese, fx];
	[tt,t,fx]=readvector(filename, 'kei');
	if ~isempty(find(~isfinite(fx))) find(~isfinite(fx)), end
	kei=[kei, fx];
end

subplot(221), 
h=plot(t,ese);
set(h(1),'linestyle','-');
set(h(2),'linestyle',':');
set(h(3),'linestyle','-.');
set(h(4),'linestyle','--');

title('6D ESE (B=0)');
set(gca,'ylim',[20.0 40.0]);
legend('periodic emode=0','periodic emode=1','insulated','free space');
subplot(223), 
h=plot(t,kei);
set(h(1),'linestyle','-');
set(h(2),'linestyle',':');
set(h(3),'linestyle','-.');
set(h(4),'linestyle','--');

title('6D KEI (B=0)');
set(gca,'ylim',[0.999 1.004]);

files=str2mat('pb1m0_3d.hdf','pb1m1_3d.hdf','ib1_3d.hdf','fb1_3d.hdf');

kei=[];
ese=[];
for i=1:size(files,1),

	filename = deblank(files(i,:));
	if exist(filename,'file')==2,
		filename=which(filename);
	else
		error(sprintf('Cannot find file %s', filename));
	end

	[tt,t,fx]=readvector(filename, 'ese');
	ese=[ese, fx];
	[tt,t,fx]=readvector(filename, 'kei');
	if ~isempty(find(~isfinite(fx))) find(~isfinite(fx)), end
	kei=[kei, fx];
end

subplot(222), 
h=plot(t,ese);
set(h(1),'linestyle','-');
set(h(2),'linestyle',':');
set(h(3),'linestyle','-.');
set(h(4),'linestyle','--');

title('6D ESE (\omega_{ci}=\omega_{pi})');
set(gca,'ylim',[20.0 40.0]);
legend('periodic emode=0','periodic emode=1','insulated','free space');
subplot(224), 
h=plot(t,kei);
set(h(1),'linestyle','-');
set(h(2),'linestyle',':');
set(h(3),'linestyle','-.');
set(h(4),'linestyle','--');

title('6D KEI (\omega_{ci}=\omega_{pi})');
set(gca,'ylim',[0.999 1.004]);

orient landscape; set(gcf,'PaperOrientation','portrait');

if exist('exportfig'),
	exportfig(gcf,'energy3d','format','eps');
else
	print -deps energy3d;
end
