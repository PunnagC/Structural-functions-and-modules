function [ structural ] = structural_nonlinearity_module(structural,panelc,seg,Np,idx_start,idx_end)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% aero.Np        = 500;
         
% keyboard

%% Error checking at INPUT
IN_min = 3; %min inputs
IN_max = 6; %max inputs
narginchk(IN_min,IN_max);
if nargin <= IN_min
    idx_start = 1;
    idx_end   = numel(panelc.rhoA);
    idx       = idx_start:1:idx_end;
else
    idx_start = 1;% sum(segc.Npanel(1:2)) + 1; %75
    idx_end   = numel(panelc.rhoA);%idx_start -1 + segc.Npanel(3); %309
    idx       = round(linspace(idx_start,idx_end,Np));
end
% keyboard
% panel_lag     = structural.lag;
% panel_bending = structural.bending;
% 
% NLbend_geom = Azrar_White_nonlinearity(panel_bending,panel,sum(seg.l));
% NLlag_geom  = Azrar_White_nonlinearity(panel_lag,panel,sum(seg.l));
% 
% structural.NLbend_geom = NLbend_geom;
% structural.NLlag_geom  = NLlag_geom;

L_tot   = sum(panelc.dr(idx_start:idx_end));
% keyboard
A       = panelc.Acs(idx);   %(idx_start:idx_end,1)
dl      = panelc.dr(idx);    %(idx_start:idx_end,1)
width   = panelc.chord(idx); %(idx_start:idx_end,1)

PHI     = structural.bending.PHI(idx,:); %(idx_start:idx_end,:)
PSI     = structural.torsion.PSI(idx,:); %(idx_start:idx_end,:)

EA      = panelc.EA(idx,1);  %(idx_start:idx_end,1)
% E      = panel.E(idx,1);  %(idx_start:idx_end,1)
EIp     = panelc.EIp(idx,1); %(idx_start:idx_end,1)
% Ip     = panel.Jp(idx,1); %(idx_start:idx_end,1) %%REMOVED PANELC to PANEL%%%%%%

dPHI_dx   = structural.bending.dPHI_dx(idx,:); %(idx_start:idx_end,:)
d2PHI_dx2 = structural.bending.d2PHI_dx2(idx,:); %(idx_start:idx_end,:)

dPSI_dx   = structural.torsion.dPSI_dx(idx,:); %(idx_start:idx_end,:)
d2PSI_dx2 = structural.torsion.d2PSI_dx2(idx,:); %(idx_start:idx_end,:)

% keyboard
cprintf('*string','*****[Drachinsky-Raveh structural nonlinearity]*****\n');
[NLh_geom,NLa_geom] = Drachinsky_Raveh_Bending_NLV3(EA,dl,width,PSI,PHI,EIp,dPHI_dx, d2PHI_dx2, dPSI_dx, d2PSI_dx2,L_tot);%,A,type,PHI
% [NLh_geom,NLa_geom] = Drachinsky_Raveh_Bending_NLV4(EA,dl,width,PSI,PHI,EIp,dPHI_dx, d2PHI_dx2, dPSI_dx, d2PSI_dx2,L_tot,panel);%,A,type,PHI
%                                                    (E,A,dl,width,PSI,PHI,Ip,dPHI_dx, d2PHI_dx2, dPSI_dx, d2PSI_dx2,L_tot)%,A,type,PHI
% [NLh_geom,NLa_geom] = Drachinsky_Raveh_nonlinearity(panelh,panela,geom,sim,aero,T,mat)%,A,type,PHI
structural.NLbend_geom     = NLh_geom;
structural.NLtorsion_geom  = NLa_geom;
cprintf('*string','*****************[Completed]*****************\n');
end

