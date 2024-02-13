function [n,u,T,S,K]=phasespace_to_moments(p)
% function [n,u,T,S,K]=phasespace_to_moments(p)

%
% $Id: phasespace_to_moments.m,v 1.4 2011/03/26 15:36:05 patrick Exp $
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

if ~isempty(findstr(p.varname,'vx')) & ...
	~isempty(findstr(p.dimsname{2},'y'))

	V=repmat(p.dims{3}',size(p.dims{2}));

	n.var=zeros(length(p.dims{1}),length(p.dims{2}));
	n.varname='Density';
	n.dims{1}=p.dims{1}; n.dims{2}=p.dims{2};
	n.dimsname{1}=p.dimsname{1}; n.dimsname{2}=p.dimsname{2};

	u.var=zeros(length(p.dims{1}),length(p.dims{2}));
	u.varname='Mean velocity';
	u.dims{1}=p.dims{1}; u.dims{2}=p.dims{2};
	u.dimsname{1}=p.dimsname{1}; u.dimsname{2}=p.dimsname{2};

	T.var=zeros(length(p.dims{1}),length(p.dims{2}));
	T.varname='Temperature';
	T.dims{1}=p.dims{1}; T.dims{2}=p.dims{2};
	T.dimsname{1}=p.dimsname{1}; T.dimsname{2}=p.dimsname{2};

	S.var=zeros(length(p.dims{1}),length(p.dims{2}));
	S.varname='Skewness';
	S.dims{1}=p.dims{1}; S.dims{2}=p.dims{2};
	S.dimsname{1}=p.dimsname{1}; S.dimsname{2}=p.dimsname{2};

	K.var=zeros(length(p.dims{1}),length(p.dims{2}));
	K.varname='Kurtosis';
	K.dims{1}=p.dims{1}; K.dims{2}=p.dims{2};
	K.dimsname{1}=p.dimsname{1}; K.dimsname{2}=p.dimsname{2};

	for t=1:length(p.dims{1})

		% density
  	n.var(t,:)=sum(squeeze(p.var(t,:,:))'); 
		ii=find(n.var(t,:)>50);

		% mean drift
		u.var(t,ii)=sum(V(:,ii).*squeeze(p.var(t,ii,:))')./n.var(t,ii);

		% temperature
		U=repmat(u.var(t,:),size(p.dims{3}'));
		T.var(t,ii)=sum((V(:,ii)-U(:,ii)).^2.*squeeze(p.var(t,ii,:))')./ ...
			n.var(t,ii);

		% skewness
		S.var(t,ii)=sum((V(:,ii)-U(:,ii)).^3.*squeeze(p.var(t,ii,:))')./ ...
			n.var(t,ii)./T.var(t,ii).^(3/2);

		% kurtosis
		K.var(t,ii)=sum((V(:,ii)-U(:,ii)).^4.*squeeze(p.var(t,ii,:))')./ ...
			n.var(t,ii)./T.var(t,ii).^2-3;

	end

end
