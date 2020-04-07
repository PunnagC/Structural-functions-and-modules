function NLM = Drachinsky_dynamic_NL_restoring_forces_moments(PHI, PSI, dz, rh, ra,h_sign, panel)%dynamic_NL(panelh, panela, rh, ra)
% h_sign is -1 when heave h is upwards (CS+)   - modified derivation
% h_sign is +1 when heave h is downwards (CS-) - usual derivation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard

x  = panel.x;
dx = panel.dx;

hp    = PHI*rh;
theta = PSI*ra;
w     = -h_sign*hp + sin(theta).*x; %Combined bending-torsion motion in physical coordinate

dz_matrix = dz*ones(1,panel.Nb);
w         = [zeros(1,panel.Nb);w]; %adding an extra row of zeros at the begining
delta_w   = diff(w,1,1); %'diff' reduces 1 row
dw_dz     = delta_w./dz_matrix; %numerical derivative

dynamic_T_term = (panel.Et.*dz)'*(dw_dz.^2); %1xnb
dynamic_T      = (0.5/panel.L)*(dynamic_T_term);
T              = 1*panel.T0_c + dynamic_T; %total tension
% keyboard
% df_NL_term = -h_sign*dynamic_T.*dw_dz;               %npxnb
df_NL_term = -h_sign*T.*dw_dz;               %npxnb
df_NL_term = [zeros(1,panel.Nb);df_NL_term]; %adding an extra row of zeros at the begining
df_NL      = diff(df_NL_term,1,1)./dz_matrix; %npxnb, %'diff' reduces 1 row
Fr_bending = ((df_NL*dx').*dz)'*PHI; %npx1
Mr_theta   = ((df_NL*(x.*dx)').*cos(theta).*dz)'*PSI; %npx1


NLM        = 1.*[Fr_bending'; Mr_theta'];

end


%%Gradient/diff function testing
% TestM = [5 3 1;
%          6 5 4;
%          7 7 7;
%          8 9 10;
%          9 11 13];
% 
%      diff(TestM,1,1)