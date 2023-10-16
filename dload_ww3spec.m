function [spec_ww3 f_ww3 dir_ww3 t_ww3 bw_ww3] = dload_ww3spec(ndbcid,varname,tlims,tres)
% -------------------------------------------------------------------------
% DLOAD_WW3SPEC  Downloads spectral model data from CDIP THREDDS 
% -------------------------------------------------------------------------
%   Syntax:
%      [spec_ww3 f_ww3 dir_ww3 t_ww3] = dload_ww3spec(ndbcid,varname,tlims,tres)
%
%   Inputs:
%      NDBCID 
%      VARNAME      options: 'EnergyDensity' or 
%                            'DirectionalSpectrum' or
%                            'Mean Direction'
%      TLIMS        default = [NaN NaN], can also use numeric indexes...
%      TRES         default = 1;
% 
%   Output:
%      [var_ww3]
% 
%   Uses:
%      get_tinfo.m (subfunction)
%
% Updated as of 08-18-2019 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
if isnumeric(ndbcid)
    ndbcid = num2str(ndbcid);
end

ww3name = ndbcid;

if ~exist('tlims')
    tlims = [NaN NaN];
elseif ~isnan(tlims)
    if tlims(2)<tlims(1)
        error('Incorrect time limits')
    end
end
if isempty(tlims)
    tlims = [NaN NaN];
end

if ~exist('tres')
    tres = 1;
end

if strcmp(varname, 'EnergyDensity') | strcmp(varname, 'MeanDirection') | contains(varname, 'Value')
    dim = 2;
elseif strcmp(varname, 'DirectionalSpectrum')
    dim = 3;
else
    error('Provided variable name is not real');
end

vo = 0;
fstr = '[0:1:49]';
fstr = '%5B0:1:49%5D'; fN = 50;
thstr = '%5B0:1:35%5D'; thN = 36;


%% LOAD WW3 from CDIP THREDDS
% -------------------------------------------------------------------------

%%%s ww3 spectral output
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/test/WW3/';
nameurl = [ww3name '_WW3_realtime.nc.ascii?'];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get frequency and directions

% Load frequency
paramurl = ['waveFrequency' fstr];
url = [baseurl nameurl paramurl]; %full url name

data = urlread(url); data = strsplit(data,'\n');
dstart = find(strcmp(data,'---------------------------------------------'))+1;

waveVar = strsplit(data{dstart+1}, ', ');
waveVar = cellfun(@str2double,waveVar);
f_ww3 = waveVar(1:end);


% Load bandwidth
paramurl = ['waveBandwidth' fstr];
url = [baseurl nameurl paramurl]; %full url name

data = urlread(url); data = strsplit(data,'\n');
dstart = find(strcmp(data,'---------------------------------------------'))+1;

waveVar = strsplit(data{dstart+1}, ', ');
waveVar = cellfun(@str2double,waveVar);
bw_ww3 = waveVar(1:end);
% bw_ww3 = NaN;


% Load direction
paramurl = ['waveDirection' thstr];
url = [baseurl nameurl paramurl]; %full url name

data = urlread(url); data = strsplit(data,'\n');
dstart = find(strcmp(data,'---------------------------------------------'))+1;

waveVar = strsplit(data{dstart+1}, ', ');
waveVar = cellfun(@str2double,waveVar);
dir_ww3 = waveVar(1:end);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be


if tlims(1)<700000
    starti = tlims(1);
    endi = tlims(2);
    Nt = endi-starti+1;
else
    [~, ~, starti, endi, Nt] = get_tinfo(baseurl, ww3name, tlims);
    Nt = endi-starti;
    Nt = Nt/tres+1;
    
end
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

if vo; disp(['     tstring for url: ' tstring]); end;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create empty variable, url name, load data, reshape

%%%% waveDirectionalEnergy %%%
if dim==3 
%     paramurl = ['wave' varname tstring '[0:1:49][0:1:35]'];
    paramurl = ['wave' varname tstring fstr thstr];
    var_ww3 = NaN(36,50*(Nt));
    url = [baseurl nameurl paramurl]; %full url name
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    for j=1:50*Nt
        waveVar = strsplit(data{dstart+j}, ', ');
        if vo; disp(waveVar{1}); end
        waveVar = cellfun(@str2double,waveVar);
        var_ww3(:,j) = waveVar(2:end);
    end
    var_ww3 = reshape(var_ww3, [36,50,Nt]);  
    
%%%% waveEnergyDensity %%%
elseif dim==2 
%     paramurl = ['wave' varname tstring '[0:1:49]'];
    paramurl = ['wave' varname tstring fstr];
    var_ww3 = NaN(50,Nt);
    url = [baseurl nameurl paramurl]; %full url name
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;
    
    for j=1:Nt
        waveVar = strsplit(data{dstart+j}, ', ');
        if vo; disp(waveVar{1}); end
        waveVar = cellfun(@str2double,waveVar);
        var_ww3(:,j) = waveVar(2:end);
    end
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
spec_ww3 = var_ww3;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% THIRD: provide time variable

tstart = find(strcmp(data,['wave' varname '.waveTime[' num2str(Nt) ']']));
waveVar = strsplit(data{tstart+1}, ', ');
waveVar = cellfun(@str2double,waveVar);
toff = datenum(1970,1,1,0,0,0);
t_ww3 = waveVar./(24*60*60) + toff;

end


function [tstart tend starti endi Nt] = get_tinfo(baseurl, ww3name, tlims)
    infourl = [baseurl ww3name '_WW3_realtime.nc.html'];
    info = urlread(infourl); info = strsplit(info,'\n');

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = Nt   
    idx = find(contains(info, 'Int32 waveTime[waveTime = ')); idx = idx(1);
    infot = strsplit(info{idx}, '= '); infot = infot{2};
    Nt = str2num(infot(1:end-2))-1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = tstart, tend
    idx = find(contains(info, 'date_created'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tstart = datenum(yr,mo,day);
    
    idx = find(contains(info, 'date_modified'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tend = datenum(yr,mo,day);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     disp(['     CDIP runs from ' datestr(tstart) ' to '  datestr(tend)])
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = starti, endi
    if ~isempty(tlims)
        endi = Nt;
        sdi = 30; % change this if you want to make the time cuttoff more accurate
        starti = Nt-floor(Nt/sdi)*sdi;
        % downloading coarse resolution time series data to get
        % approximate index location for tlims
        endi = floor(endi/sdi).*sdi;
        tstring = [ '[' num2str(starti) ':' num2str(sdi) ':' num2str(endi) ']']; 
        tstring = [ '%5B' num2str(starti) ':' num2str(sdi) ':' num2str(endi) '%5D']; 
        idx = [1 starti:sdi:endi];
        
        nameurl = [ww3name '_WW3_realtime.nc.ascii?'];
        paramurl = ['waveTime' tstring]; % to do figure out subset!
        url = [baseurl nameurl paramurl];
        waveVar = dload(url);
        
        tstring = '[1]';
        tstring = '%5B1%5D';
        paramurl = ['waveTime' tstring]; % to do figure out subset!
        url = [baseurl nameurl paramurl];
        initwaveVar = dload(url);
        
        waveVar = [initwaveVar waveVar];
        
        toff = datenum(1970,1,1,0,0,0);
        simpletime = waveVar./(24*60*60) + toff; % the course resolution data, corresponding indexes in VAR=idx

%         [~, idxstart] = min(abs(simpletime-tlims(1)));
%         [~, idxend] = min(abs(simpletime-tlims(2)));

        [~, ss] = sort(abs(simpletime-tlims(1)));
        idxstart = min(ss(1:2));
        
        [~, ss] = sort(abs(simpletime-tlims(2)));
        idxend = max(ss(1:2));

        starti = idx(idxstart); endi = idx(idxend);
    else % get the entire thing
        starti = 1;
        endi = Nt;
    end
end


function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end
