function H=channel2tf(h,delays,freqs);

% H=channel2tf(h,delays,freqs);
%
% Convert a channel specified as an impulse response into a frequency
% response measured at the frequencies specified in 'freqs'.
%
% The expected format of 'h' is n_ue_ant*n_bs_ant*n_delays.
% The format of H is n_ue_ant*n_bs_ant*n_freqs

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
error(chk_param(h,'H','numeric'));
error(chk_param(delays,'delays','real','vector','horizontal'));
error(chk_param(freqs,'freqs','vector','real'));

n_ant_ue=size(h,1);
n_ant_bs=size(h,2);

H=zeros(n_ant_ue,n_ant_bs,length(freqs));
for k=1:length(delays)
  H(:,:,:)=H(:,:,:)+reshape(kron(exp(-j*2*pi*freqs*delays(k)),h(:,:,k)),n_ant_ue,n_ant_bs,length(freqs));
end

