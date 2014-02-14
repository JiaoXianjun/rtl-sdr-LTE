function y=cyc_filt(x,h);

% y=cyc_filt(x,h)
%
% Cyclically filter signal x with filter h. h can be either an impulse
% response or a function which takes as input a value between -1 and 1
% and returns the response of the filter at that frequency.

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

error(nargchk(2,2,nargin));
error(chk_param(x,'x','vector'));
if (isa(h,'function_handle'))
  f=linspace(0,2,length(x)+1);
  f=f(1:end-1);
  y=ifft(fft(x).*h(mod(f+1,2)-1));
  return

elseif (isvector(h))
  fft_size=max(length(x),length(h));
  y=ifft(fft(x,fft_size).*fft(h,fft_size));
  return

else
  error('Unrecognizable transfer function specified');

end

% Should never reach here...
error('Check code...');

