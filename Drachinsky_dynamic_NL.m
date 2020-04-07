function NLM = Drachinsky_dynamic_NL(PHI, PSI, dPHI_dx, d2PHI_dx2, dPSI_dx, d2PSI_dx2, t0 , dx, int, rh, ra,h_sign)%dynamic_NL(panelh, panela, rh, ra)
% flag is +1 when heave h is upwards
% flag is -1 when heave h is downwards
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h_x  = -h_sign*dPHI_dx*rh;      %Npx1 size
% h_xx = -h_sign*d2PHI_dx2*rh;    %Npx1 size 

h_x  = -h_sign*dPHI_dx*rh;      %Npx1 size
h_xx = -h_sign*d2PHI_dx2*rh;    %Npx1 size 

theta_x  = dPSI_dx*ra;   %Npx1 size
theta_xx = d2PSI_dx2*ra; %Npx1 size


%% Evaluating 1st part of BDOF tB1

int_dy    = int.int_dy;
int_y_dy  = int.int_y_dy;
int_y2_dy = int.int_y2_dy;
int_y3_dy = int.int_y3_dy;
int_y4_dy = int.int_y4_dy;


intx_tB11    = 2*dx'*(t0.*h_x.*h_xx);
intx_tB12    = 2*dx'*(t0.*(h_x.*theta_xx + theta_x.*h_xx));
intx_tB13    = 2*dx'*(t0.*theta_x.*theta_xx);

inty_tB11_1 = (intx_tB11.*h_x);
inty_tB11_2 = (intx_tB11.*theta_x);
inty_tB11   = inty_tB11_1*int_dy + inty_tB11_2*int_y_dy;

inty_tB12_1 = (intx_tB12.*h_x);
inty_tB12_2 = (intx_tB12.*theta_x);
inty_tB12   = inty_tB12_1*int_y_dy + inty_tB12_2*int_y2_dy;

inty_tB13_1 = (intx_tB13.*h_x);
inty_tB13_2 = (intx_tB13.*theta_x);
inty_tB13   = inty_tB13_1*int_y2_dy + inty_tB13_2*int_y3_dy;

int_tB1  = (inty_tB11 + inty_tB12 + inty_tB13);
% keyboard
%% Evaluating 2nd part of BDOF tB2
intx_tB21 = dx'*(t0.*h_x.^2);
intx_tB22 = dx'*(t0.*.2.*h_x.*theta_x);
intx_tB23 = dx'*(t0.*theta_x.^2);

int_tB2 = (intx_tB21*int_dy + intx_tB22*int_y_dy + intx_tB23*int_y2_dy)*h_xx ...
          + (intx_tB21*int_y_dy + intx_tB22*int_y2_dy + intx_tB23*int_y3_dy)*theta_xx;



NLb       = ((int_tB1 + int_tB2).*dx)'*PHI; %for bending

%% Evaluating the torsion part

inty_tT11 = inty_tB11_1*int_y_dy + inty_tB11_2*int_y2_dy;
inty_tT12 = inty_tB12_1*int_y2_dy + inty_tB12_2*int_y3_dy;
inty_tT13 = inty_tB13_1*int_y3_dy + inty_tB13_2*int_y4_dy;

int_tT1  = (inty_tB11 + inty_tB12 + inty_tB13);

int_tT2 = (intx_tB21*int_y_dy + intx_tB22*int_y2_dy + intx_tB23*int_y3_dy)*h_xx ...
          + (intx_tB21*int_y2_dy + intx_tB22*int_y3_dy + intx_tB23*int_y4_dy)*theta_xx;

NLa       = ((int_tT1 + int_tT2).*dx)'*PSI; %for torsion
% keyboard
NLM = [NLb.*1, NLa.*1]';
end