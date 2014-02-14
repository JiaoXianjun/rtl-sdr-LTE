function tfplot(imp);

% tfplot(imp);
%
% Plot the transfer function of the filter whose impulse response is imp.
%
% Multiple transfer functions can be plotted simultaneously where each
% impulse response represents one row of imp.

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

error(nargchk(1,1,nargin));
%error(chk_param(imp,'imp','vector'));

n_pts_min=4096;
if (isvector(imp))
  n_pts=max(n_pts_min,length(imp));
else
  n_pts=max(n_pts_min,size(imp,2));
end

f=linspace(0,2,n_pts+1);
f=f(1:end-1);
f(f>=1)=f(f>=1)-2;
f=fftshift(f);
m=transpose(fftshift(fft(transpose(imp),n_pts)));
plot(f,db20(abs(m)));

xlabel('Frequency/(fs/2)');
ylabel('Power gain (dB)');
title('Transfer function');
zgo;

