function state=etc(a1,a2)

% Simple function to estimate time to completion (ETC).
%
% Usage:
%
% state=etc;
% for t=1:n_trials
%   do_something(t);
%   state=etc(state,t/n_trials);
% end

% Copyright 2012 Evrytania LLC (http://www.evrytania.com)
%
% Written by James Peroulas <james@evrytania.com>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Specified in seconds
update_delay=5;

error(nargchk(0,2,nargin));
if (nargin==1)
  error('Either 0 or 2 arguments can be supplied to this function');
end
if (nargin==0)
  state.sim_start=clock;
  state.last_print=clock;
  return
end

error(chk_param(a1,'state','struct'));
state=a1;
error(chk_param(a2,'completion','scalar','real','>',0,'<=',1));
completion=a2;

if (completion==1)
  ttime=round(etime(clock,state.sim_start));
  disp(['Total elapsed time: ' time_print(ttime)]);
elseif (etime(clock,state.last_print)>5)
  ttime=etime(clock,state.sim_start);
  remaining=round(ttime/completion-ttime);
  ttime=round(ttime);
  disp([sprintf('%4.1f%%',completion*100) ' |Elapsed: ' time_print(ttime) '| |Remaining: ' time_print(remaining) '|']);
  state.last_print=clock;
end

function str=time_print(e);

h=floor(e/3600);
e=e-h*3600;
m=floor(e/60);
e=e-m*60;
s=e;

if (h)
  str=sprintf('%dh %2dm %2ds',h,m,s);
elseif (m)
  str=sprintf('%2dm %2ds',m,s);
else
  str=sprintf('%2ds',s);
end

