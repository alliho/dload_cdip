function [lat, lon] = dload_cdiplatlon(cdipid)
% -------------------------------------------------------------------------
% dload_cdiplatlon  Downloads model data from CDIP THREDDS (historic and rt)
% -------------------------------------------------------------------------
%   Syntax:
%      [lat, lon] = dload_cdiplatlon(cdipid)
%
%   Inputs:
%      CDIPID 
% 
%   Output:
%      [lat, lon]
% 
% Updated as of 08-20-2019 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 

if isnumeric(cdipid)
    cdipid = num2str(cdipid);
    if length(cdipid)<3
        cdipid = ['0' cdipid];
    end
end
cdipname = cdipid;

if ~exist('tlims')
    tlims = [];
end

if ~exist('tres')
    tres = 1;
end

%% load
try
    baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
    nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.html'];
    % http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/191p1/191p1_historic.nc.html
    url = [baseurl nameurl];
    data = urlread(url); data = strsplit(data,'\n');
catch
    baseurl = 'https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/';
    nameurl = [cdipname 'p1_rt.nc.html'];
    url = [baseurl nameurl];
    data = urlread(url); data = strsplit(data,'\n');
end

idx = find(contains(data, 'geospatial_lat_max')); 
var = strsplit(data{idx}, ': ');
lat = str2num(var{2});

idx = find(contains(data, 'geospatial_lon_max')); 
var = strsplit(data{idx}, ': ');
lon = str2num(var{2});


end