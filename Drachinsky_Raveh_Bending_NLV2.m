function  [NLh_geom,NLa_geom] = Drachinsky_Raveh_Bending_NLV2(EA,dl,width,PSI,PHI,EIp,dPHI_dx, d2PHI_dx2, dPSI_dx, d2PSI_dx2,L_tot)%,A,type,PHI
%                                                            (EA,dl,width,PSI,PHI,EJp,dPHI_dl, d2PHI_dl2, dPSI_dl, d2PSI_dl2,L_total)
%This calculates the NL matrices as per [Drachinsky Raveh]Limit-cycle oscillations of a pre-tensed membrane strip
%   Detailed explanation goes here   panelh,geom,sim,aero,T,mat,Acs,type,nmodes,PHI

[~,nh] = size(PHI);
[~,na] = size(PSI);

rh  = sym('rh', [1,nh]); % bending DOF
ra  = sym('ra',[1,na]);% torsion DOF

w_x_sym   = dPHI_dx*transpose(rh);%transpose
w_xx_sym  = d2PHI_dx2*transpose(rh);

th_x_sym  = dPSI_dx*transpose(ra);%transpose
th_xx_sym = d2PSI_dx2*transpose(ra);

Iww       = sum((w_x_sym.^2).*dl);%sum
Ithth     = sum((th_x_sym.^2).*dl);
Iwth      = sum((w_x_sym.*th_x_sym).*dl);

% For bending

% keyboard
t0       = 1*0.5/L_tot; %t0       = 0.5/L_tot;
t1       = transpose(t0.*(EA.*Iww + EIp.*Ithth).*w_xx_sym.*dl);%transpose
t2       = transpose((2*t0.*EIp.*Iwth).*th_xx_sym.*dl);%transpose
NLh_geom = (t1 + t2)*PHI;

NLh_geom = vpa(NLh_geom,3);
NLh_geom = simplify(expand(NLh_geom));
NLh_geom = matlabFunction(-NLh_geom); % NOTE the -ve sign
cprintf('text','Bending Non-Linear structural terms evaluated.\n')

% For torsion
ta1       = transpose(t0.*(EIp.*Iww + EA.*width.^4.*Ithth./80).*th_xx_sym.*dl);
ta2       = transpose((2*t0.*EIp.*Iwth).*w_xx_sym.*dl);
NLa_geom  = (ta1 + ta2)*PSI;

NLa_geom  = vpa(NLa_geom,3);
NLa_geom  = simplify(expand(NLa_geom));
NLa_geom  = matlabFunction(-NLa_geom); % NOTE the -ve sign
cprintf('text','Torsional Non-Linear structural terms evaluated.\n')
% disp('All Drachinsky-Raveh Non-Linear structural terms evaluated.')
end

