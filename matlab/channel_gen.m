function [h h_timestamp delays]=channel_gen(sim_length,n_ant_bs,n_ant_ue,v_mobile,fc,pos_bs,pos_ue,scenario);

% [h h_timestamp delays]=channel_gen(sim_length,n_ant_bs,n_ant_ue,v_mobile,fc,pos_ue,pos_bs,scenario);
%
% Produces channel responses between UE's and basestations. Either the number
% of UE's must be 1 or the number of basestations must be 1.
%
% sim_length will cause enough channel samples to be produced to simulate
%   sim_length seconds worth of time.
% n_ant_bs is the number of antennas in the basestation
% n_ant_ue is the number of antennas in the UE
% v_mobile is the velocity and direction of travel of the UE's
% fc is the carrier frequency
% pos_bs is the location of the basestations
% pos_ue is the location of the UE's
%
% scenario can be one of the following:
%   IMT_Advanced channel models:
%     InH' 'UMi' 'SMa' 'UMa' 'UMiO2I' 'RMa'
%   Basic channel models:
%     'flat' 'rayleigh'
%
% n_links channels will be created and each channel will have n_delays
% tap delays. n_delays varies based on the channel model.
% n_chan_samp samples of the channel will be created.
%
% delays is of size:
%   n_links*n_delays
% h_timestamp is of size:
%   1*n_chan_samp
% h is of size:
%   n_ant_ue*n_ant_bs*n_delays*n_chan_samp*n_links
%
% h_timestamp(3) indicates that h(:,:,:,3,:) should be used starting at
% simulation time h_timestamp(3).

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

% Constants
% Number of channels produced per half-wavelength.
sample_density=10;

error(nargchk(8,8,nargin));
error(chk_param(sim_length,'sim_length','scalar','real','>=',0));
error(chk_param(n_ant_bs,'n_ant_bs','scalar','real','integer','>=',1));
error(chk_param(n_ant_ue,'n_ant_ue','scalar','real','integer','>=',1));
error(chk_param(v_mobile,'v_mobile','vector'));
error(chk_param(fc,'fc','scalar','real','>',0));
error(chk_param(pos_ue,'pos_ue','vector','integer'));
error(chk_param(pos_bs,'pos_bs','vector','integer'));
error(chk_param(scenario,'scenario','string'));
if (length(v_mobile)~=length(pos_ue))
  error('v_mobile must be same size as pos_ue');
end

% Derive some parameters
n_bs=length(pos_bs);
n_ue=length(pos_ue);
n_links=max(n_bs,n_ue);

% Secondary checking
if ((n_bs~=1)&&(n_ue~=1))
  error('Either number of BS OR number of MS must be equal to 1');
end

% Determine how many channel samples will be created.
wavelength=speed_of_light/fc;
if (max(abs(v_mobile))==0)
  num_time_samples=1;
  h_timestamp=0;
else
  num_time_samples=ceil(max(abs(v_mobile))*sim_length/(wavelength/2)*sample_density+1);
  h_timestamp=(0:num_time_samples-1)*(wavelength/2)/sample_density/max(abs(v_mobile));
end

% Determine which channel model group is requested.
IMTA_models={'InH','UMi','SMa','UMa','UMiO2I','RMa'};
basic_models={'flat','Rayleigh','EPA','EVA','ETU'};
if (any(strcmpi(scenario,IMTA_models)))
  cm_group='IMT-A';
  scenario=find(strcmpi(scenario,IMTA_models));
elseif (any(strcmpi(scenario,basic_models)))
  cm_group='basic';
else
  error('Unrecognized channel model requested.');
end

if (strcmpi(cm_group,'basic'))
  if (strcmpi(scenario,'flat')||strcmpi(scenario,'rayleigh'))
    delays=0;
    power=1;
  elseif strcmpi(scenario,'EPA')
    delays=[0 30 70 90 110 190 410]*1e-9;
    power=udb10([0.0 -1.0 -2.0 -3.0 -8.0 -17.2 -20.8]);
  elseif (strcmpi(scenario,'EVA'))
    delays=[0 30 150 310 370 710 1090 1730 2510]*1e-9;
    power=udb10([0.0 -1.5 -1.4 -3.6 -0.6 -9.1 -7.0 -12.0 -16.9]);
  elseif (strcmpi(scenario,'ETU'))
    delays=[0 50 120 200 230 500 1600 2300 5000]*1e-9;
    power=udb10([-1.0 -1.0 -1.0 0.0 0.0 0.0 -3.0 -5.0 -7.0]);
  else
    error('Check code...');
  end

  % Some post-processing
  power=power/sum(power);
  n_delays=length(delays);
  delays=repmat(delays,n_links,1);

  % Build the response
  h=NaN(n_ant_ue,n_ant_bs,n_delays,num_time_samples,n_links);
  for t=1:n_links
    if (n_ue==1)
      v=v_mobile;
    else
      v=v_mobile(t);
    end
    dop_shift=doppler_shift(fc,v);
    for k=1:n_delays
      for l=1:n_ant_ue
        for m=1:n_ant_bs
          if (strcmpi(scenario,'flat'))
            taps=ones(1,num_time_samples);
          else
            taps=sqrt(power(k))*rayleigh(num_time_samples,2*dop_shift,1/((wavelength/2)/sample_density/v)+eps);
          end
          h(l,m,k,:,t)=taps;
        end
      end
    end
  end

elseif (strcmpi(cm_group,'IMT-A'));
  % Select the orientation of the antenna arrays
  if (n_bs==1)
    bs_array_ang=0;
    ue_array_ang=(2*rand(1,n_ue)-.5)*pi;
  else
    ue_array_ang=0;
    bs_array_ang=(2*rand(1,n_bs)-.5)*pi;
  end
  % Angle between the basestation and the ue
  bs_ue_ang=angle(pos_ue-pos_bs);
  ue_bs_ang=angle(pos_bs-pos_ue);

  % There's a bug in wim if n_ant_ue or n_ant_bs is 1... Work around this...
  n_ant_ue_orig=n_ant_ue;
  n_ant_bs_orig=n_ant_bs;
  n_ant_ue=max([2 n_ant_ue]);
  n_ant_bs=max([2 n_ant_bs]);

  % Another bug is if v_mobile is zero
  v_mobile(v_mobile==0)=.000001;

  % wimpar
  wimpar=struct( ...
    'NumBsElements',n_ant_bs,...
    'NumMsElements',n_ant_ue,...
    'range',1,...
    'end_time',1,... % Observation end time for B5 - time points are taken as:  wimpar.TimeVector=linspace(0,wimpar.end_time,T);
    'SampleDensity', sample_density,...   % in samples/half-wavelength
    'NumTimeSamples',num_time_samples,...
    'UniformTimeSampling','yes',...
    'IntraClusterDsUsed','yes',...        % Two strongest clusters are divided into three subclusters
    'NumSubPathsPerPath',20,...           % only value supported is 20.
    'FixedPdpUsed','no',...               % Use fixed delays and path powers
    'FixedAnglesUsed','no',...            % Use fixed AoD/AoAs
    'PolarisedArrays','no',...           % use polarised arrays
    'TimeEvolution','no',...              % use of time evolution option
    'CenterFrequency',fc,...              % in Herz
    'DelaySamplingInterval',0,...
    'PathLossModelUsed','no',...
    'ShadowingModelUsed','no',...
    'PathLossModel','pathloss',...
    'AnsiC_core','no',...
    'LookUpTable',0,...                   % number of points in Ansi-C core look-up table for cosine, 0 if not used
    'RandomSeed',[],...                   % if empty, seed is not set.
    'UseManualPropCondition','yes' ...
  );
  % linkpar
  linkpar=struct( ...
    'ScenarioVector',repmat(scenario,1,n_links),... % A2, B1, B4, C1, C2, D1
    'PropagConditionVector',zeros(1,n_links),...
    'MsBsDistance',abs(pos_ue-pos_bs),...
    'BsHeight',repmat(NaN,1,n_links),... % NaN -> default heights from D1.1.2
    'MsHeight',repmat(NaN,1,n_links),... % NaN -> default heights from D1.1.2
    'ThetaBs',180/pi*(bs_array_ang-bs_ue_ang),...
    'ThetaMs',180/pi*(ue_array_ang-ue_bs_ang),...
    'MsVelocity',abs(v_mobile),...
    'MsDirection',180/pi*(pi/2-angle(v_mobile)),...
    'StreetWidth',20*ones(1,n_links),...
    'LayoutType',0*ones(1,n_links),...  % Layout type for UMi (B1/C4) path loss, 0=hexagonal, 1=Manhattan
    'NumFloors',1*ones(1,n_links),... The ground floor is number 1
    'OtoI_OutdoorPL',1*ones(1,n_links),... Outdoor-to-Indoor propagation condition for UMi, 1 is LOS, 0 is NLOS (see. note 3 in A1-2, [4])
    'Dist1',repmat(NaN,1,n_links),...
    'BuildingHeight', repmat(NaN,1,n_links),... % NaN -> default heights from ScenParTables
    'LoS02ILinks', repmat(0,1,n_links),... % NaN -> will be drawn randomly
    'LoS02VLinks', repmat(0,1,n_links),... % NaN -> will be drawn randomly
    'NLoS02ILinks', repmat(0,1,n_links),... % NaN -> will be drawn randomly
    'NLoS02VLinks', repmat(0,1,n_links) ... % NaN -> will be drawn randomly
  );
  % antpar
  antpar=struct( ...
    'BsGainPattern',{1},...                         % in general: [Number_of_antennas, 2, Elevation_points, Azimuth_points]
    'BsGainAnglesAz',{linspace(-180,176,90)},...    % size [1 Azimuth_points]
    'BSGainAnglesEl',{0},...                        % size [1 Elevation_points] (parameter ignored)
    'BsElementPosition',[0.5],...                   % in wavelengths. When scalar, uniform spacing assumed
    'MsGainPattern',{1},...
    'MsGainAnglesAz',{linspace(-180,176,90)},...
    'MsGainAnglesEl',{0},...
    'MsElementPosition',[0.5],...                   % in wavelengths. When scalar, uniform spacing assumed
    'InterpFunction','interp_gain',...              % name of the interpolation function
    'InterpMethod','cubic' ...                    % interpolation method, depends on the function used
  );

  % Generate the channels!
  [h delays]=wim(wimpar,linkpar,antpar);
  % Fix the bug in wim...
  h=h(1:n_ant_ue_orig,1:n_ant_bs_orig,:,:,:);

else
  error('Check code...');

end

