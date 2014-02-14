function ch=rayleigh(n_samp,fd,fs)

% ch=rayleigh(n_samp,fd,fs)
%
% Create n_samp samples of a Rayleigh channel whose doppler spread is fd
% and which is sampled at frequency fs.

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

error(nargchk(3,3,nargin));
error(chk_param(n_samp,'n_samp','scalar','real','integer','>=',1));
error(chk_param(fd,'fd','scalar','real','>=',0));
error(chk_param(fs,'fs','scalar','real','>',0));
if (fd>fs/2)
  error('Doppler spread must be smaller than fs/2');
end

% Special case when requested doppler spread is zero.
if (fd==0)
  ch=rayleigh(1,0.5,2);
  ch=repmat(ch,1,n_samp);
  return;
end

% First, a rayleigh channel is created with a doppler frequency of 0.5
% and a sampling frequency of 2. Then, this channel is interpolated
% up to the desired sampling frequency.

% Calculate the interpolation factor.
intp=fs/fd/4;

% Rayleigh transfer function
rtf=@(f)(abs(f)<.5).*sqrt(1./((pi*0.5)*sqrt(1-(f/(0.5+eps(0.5))).^2)));

% Create the source Rayleigh signal
n_samp_source=max(2048,ceil(n_samp/intp));
bb=sqrt(2)*cyc_filt(blnoise(n_samp_source),rtf);

% Now, interpolate up to the desired sampling rate.
ch=interpft(bb,round(n_samp_source*intp));
ch=ch(1:n_samp);

