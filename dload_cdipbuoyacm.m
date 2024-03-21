function [cdip] = dload_cdipbuoyacm(cdipid,tlims, tres)
% -------------------------------------------------------------------------
% DLOAD_CDIPBUOYACM  Downloads all info about CDIP buoy
% -------------------------------------------------------------------------
%   Syntax:
%      [cdip] = dload_cdipbuoyacm(cdipid,tlims, tres)
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
% 
% Updated as of 03-12-2021 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


varname = 'acmTime';
cdip.time  = dload_cdipvar(cdipid, varname, tlims, tres);
varname = 'acmSpeed';
cdip.speed  = dload_cdipvar(cdipid, varname, tlims, tres);
varname = 'acmDirection';
cdip.dir  = dload_cdipvar(cdipid, varname, tlims, tres);

r = cdip.speed; theta = deg2rad(cdip.dir);
cdip.u = r.*cos(theta); cdip.v = r.*sin(theta);



