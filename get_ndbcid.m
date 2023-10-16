function [ndbcid] = get_ndbcid(cdipid)
% -------------------------------------------------------------------------
% get_ndbcid  Gets NDBC ID from CDIP ID
% -------------------------------------------------------------------------
%   Syntax:
%      [ndbcid] = get_ndbcid(cdipid)
%
%   Inputs:
%      CDIPID 
% 
%   Output:
%      [ndbcid]
% 
% Updated as of 01-23-2020 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
idstr = cdipid;

if isnumeric(idstr)
    idstr = num2str(idstr);
    if length(idstr)<3
        idstr = ['0' idstr];
    end
end


%% load
try
    url = ['http://cdip.ucsd.edu/m/products/?stn=' idstr 'p1&param=waveHs'];
    dat = urlread(url); dat = strsplit(dat,'\n');
    idx = find(contains(dat, 'NDBC')); 
    % tmp = strsplit(dat{idx+1}, ',');
    tmp = strsplit(dat{idx}, ':');
    tmp = regexp(tmp{2}, 'href=(\S+)(\s*)$', 'tokens', 'lineAnchors');

    tmp = strsplit(tmp{1}{1}, '=');
    tmp = tmp{2};
    ndbcid = str2num(tmp(1:5));

catch
    url = ['http://cdip.ucsd.edu/themes/cdip/?pb=1&u2=s:' idstr ':st:1&d2=p9'];
    dat = urlread(url); dat = strsplit(dat,'\n');
    idx = find(contains(dat, 'NDBC')); 
    tmp = strsplit(dat{idx}, 'NDBC');
    tmp = strsplit(tmp{2}, '<');
    tmp = tmp{1};
    ndbcid = str2num(tmp);
end




end
