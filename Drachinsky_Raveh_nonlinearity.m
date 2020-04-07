function  [NLh_geom,NLa_geom] = Drachinsky_Raveh_nonlinearity(panelh,panela,geom,sim,aero,T,mat)%,A,type,PHI
%This calculates the NL matrices as per [Drachinsky Raveh]Limit-cycle oscillations of a pre-tensed membrane strip
%   Detailed explanation goes here   panelh,geom,sim,aero,T,mat,Acs,type,nmodes,PHI
[~,nh] = size(panelh.PHI);
[~,na] = size(panela.PSI); %finding the number of torsional modes

Acs   = panelh.Acs;
dl    = panelh.dl;
width = panelh.width_vec;
PSI   = panela.PSI;
PHI   = panelh.PHI;

syms x
q  = sym('q', [1,nh]); % bending DOF
r  = sym('r',[1,na]);% torsion DOF
Ip = panela.gamma_vec;
% Ip = panela.J_vec;

w_x_sym   = panelh.dPHI_dl*transpose(q);%transpose
w_xx_sym  = panelh.d2PHI_dl2*transpose(q);

th_x_sym  = panela.dPSI_dl*transpose(r);%transpose
th_xx_sym = panela.d2PSI_dl2*transpose(r);

Iww       = sum((w_x_sym.^2).*dl);%sum
Ithth     = sum((th_x_sym.^2).*dl);
Iwth      = sum((w_x_sym.*th_x_sym).*dl);

% For bending
t0       = 3*0.5*(panelh.E_vec)/sum(geom.L); %NOTE: PRE-MULTIPLIED with 3
t1       = transpose(t0.*(Acs.*Iww + Ip.*Ithth).*w_xx_sym.*dl);
t2       = transpose((2*t0.*Ip.*Iwth).*th_xx_sym.*dl);
NLh_geom = (t1 + t2)*PHI;

NLh_geom = vpa(NLh_geom,3);
NLh_geom = simplify(expand(NLh_geom));
NLh_geom = matlabFunction(-NLh_geom); % NOTE the -ve sign

% For torsion
ta1       = transpose(t0.*(Ip.*Iww + Acs.*width.^4.*Ithth./80).*th_xx_sym.*dl);
ta2       = transpose((2*t0.*Ip.*Iwth).*w_xx_sym.*dl);
NLa_geom  = (ta1 + ta2)*PSI;

NLa_geom  = vpa(NLa_geom,3);
NLa_geom  = simplify(expand(NLa_geom));
NLa_geom  = matlabFunction(-NLa_geom); % NOTE the -ve sign

end

