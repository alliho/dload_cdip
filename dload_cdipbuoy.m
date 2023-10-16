function [cdip] = dload_cdipbuoy(cdipid,tlims, tres, varargin)
% -------------------------------------------------------------------------
% DLOAD_CDIPVAR  Downloads all info about CDIP buoy
% -------------------------------------------------------------------------
%   Syntax:
%      [cdip] = dload_cdipbuoy(cdipid,tlims, tres)
%
%   Inputs:
%      CDIPID 
%      TLIMS        default = [Nt-100 to Nt];
%      TRES         default = 1;
% 
%   Output:
%      [cdip]
% 
%   Uses:
%      get_tinfo.m (subfunction)
%      a lot of other things
% 
%   Sample: 
%      [cdip] = dload_cdipbuoy(cdipid, [t1 t2],1, 'include', [{'hs' 'tp' 'dp' 'md'}])
% 
% Updated as of 03-12-2021 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% OPTIONS
%%% varargin - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
omiopt = 0;
incopt = 0;
for i=1:length(varargin)
  vin = varargin{i};
  if isequal(vin,'omit')
    opt = varargin{i+1};
    omi = opt;
    if ~isempty(omi)
        omiopt = 1;
    end
  end
  if isequal(vin,'include')
    opt = varargin{i+1};
    inc = opt;
    incopt = 1;
  end
end

%% LOAD
cdip.depth = dload_cdipdepth(cdipid);
[cdip.lat0, cdip.lon0] = dload_cdiplatlon(cdipid);
% -------------------------------------------------------------------------
varname = 'Time';
cdip.time  = dload_cdipvar(cdipid, varname, tlims, tres);



varnms = [{'Hs', 'Tp', 'gpsLatitude', 'gpsLongitude', 'Dp'}];
fldnms = [{'hs', 'tp', 'lat', 'lon', 'dp'}];
fldnms = [{'hs', 'tp', 'lat', 'lon', 'dp'}];

if omiopt
    idx = find(~contains(fldnms,omi));
    varnms = varnms(idx);
    fldnms = fldnms(idx);
elseif incopt
    idx = find(contains(fldnms,inc));
    varnms = varnms(idx);
    fldnms = fldnms(idx);
end

for i=1:length(varnms)
    varname = varnms{i};
    fldnm = fldnms{i};
    disp(['Downloading: ' varname])
    cdip.(fldnm) = dload_cdipvar(cdipid, varname, tlims, tres);
end
% varname = 'Time';
% cdip.time  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'Hs';
% cdip.hs  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'Tp';
% cdip.tp  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'gpsLatitude';
% cdip.lat  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'gpsLongitude';
% cdip.lon  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'Dp';
% cdip.dp  = dload_cdipvar(cdipid, varname, tlims, tres);

% varname = 'acmTime';
% cdip.current.time  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'acmSpeed';
% cdip.current.speed  = dload_cdipvar(cdipid, varname, tlims, tres);
% varname = 'acmDirection';
% cdip.current.dir  = dload_cdipvar(cdipid, varname, tlims, tres);
% 
% r = cdip.current.speed; theta = deg2rad(cdip.current.dir);
% cdip.current.u = r.*cos(theta); cdip.current.v = r.*sin(theta);
% 



varnms = [{'EnergyDensity', 'MeanDirection', 'CheckFactor', 'A1Value', 'B1Value', 'A2Value', 'B2Value'}];
fldnms = [{'sf', 'md', 'check', 'a1', 'b1', 'a2', 'b2'}];


if omiopt
    idx = find(~contains(fldnms,omi));
    varnms = varnms(idx);
    fldnms = fldnms(idx);
elseif incopt
    idx = find(contains(fldnms,inc));
    varnms = varnms(idx);
    fldnms = fldnms(idx);
end



for i=1:length(varnms)
    varname = varnms{i};
    fldnm = fldnms{i};
    disp(['Downloading: ' varname])
    [cdip.(fldnm) cdip.f cdip.timesf] = dload_cdipspec(cdipid,varname,tlims,tres);
    cdip.df = gradient(cdip.f);

end


% varname = 'EnergyDensity';
% [cdip.sf cdip.f ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% varname = 'MeanDirection';
% [cdip.md, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% 
% %%% additional shit
% varname = 'CheckFactor';
% [cdip.check, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% 
% varname = 'A1Value';
% [cdip.a1, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% 
% varname = 'B1Value';
% [cdip.b1, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% 
% varname = 'A2Value';
% [cdip.a2, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);
% 
% varname = 'B2Value';
% [cdip.b2, ~, ~] = dload_cdipspec(cdipid,varname,tlims,tres);


% cdip.sigma1 = sqrt(2.*(1-sqrt(cdip.a1.^2 + cdip.b1.^2)));

