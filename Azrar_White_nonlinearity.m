function [ NL_geom ] = Azrar_White_nonlinearity(struc,panel,L)
%This calculates the NL bij matrix as given in Azrar White journal paper
 

E  = panel.E;
A  = panel.Acs;
dr = panel.dr; %strip/panel width

PHI       = struc.PHI;
dPHI_dx   = struc.dPHI_dx;
d2PHI_dx2 = struc.d2PHI_dx2;

[~,nmode] = size(PHI);
q_vec     = sym('q', [1,nmode]);
disp('Here1')
w_x_sym  = dPHI_dx*transpose(q_vec);
w_xx_sym = d2PHI_dx2*transpose(q_vec);

term1    = 0.5*(E.*A)/L;
term2    = sum((w_x_sym.^2).*dr);
term3    = w_xx_sym*term2.*dr;
term3    = vpa(term3,5);
disp('Here2')

NL_geom    = PHI'*(term3.*term1);%
NL_geom    = vpa(NL_geom,3);
NL_geom    = simplify(expand(NL_geom));
NL_geom    = matlabFunction(NL_geom);
NL31_t1    = vpa(term3.*term1,3);
NL31_t2    = simplify(expand(NL31_t1));
panel.NL31 = matlabFunction(NL31_t2./dr);
% keyboard
end

