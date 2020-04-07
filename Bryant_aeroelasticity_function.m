% function data_out =  Bryant_aeroelasticity_function(T_multiplier, L_s, thk_s, AR, PZT_l_ratio, PZT_t_ratio,U_speed,a,e,h_sign)
clc;
clearvars;
close all;
% restoredefaultpath
addpath('G:\My Drive\Codes\TMM-structure-functions\');
addpath('G:\My Drive\Codes\Postprocessing functions\');
addpath('G:\My Drive\Codes\Aerodynamics-functions\Linear\Peters\');
addpath('G:\My Drive\Codes\Aerodynamics-functions\Nonlinear\Bryant model\'); 
addpath('G:\My Drive\Codes\Experimental results for code\');
addpath('G:\My Drive\Codes\misc-functions\');
addpath('G:\My Drive\Codes\Panel discretization functions\')
addpath('G:\My Drive\Codes\Structural functions and modules\');
addpath('G:\My Drive\Punnag journal paper work\[5][Progress]TensiononPiezoCoupling\Matlab code\structural parameter files\');


L_s          = 235e-3;%350e-3; %ribbon length 236e-3;%
thk_s        = 0.2e-3; %ribbon thickness 0.2e-3;%
T_multiplier = 0 ; %0.05 0.1 0.2
AR           = [10];  %aspect ratio 7  9.291;%
PZT_l_ratio  = [41/(L_s*1e3)]; %1 side pzt length to overall ribbon length ratio 0.1737;%
PZT_t_ratio  = [thk_s*1e3/0.127]; %pzt thickness to ribbon thickness ratio 1.5748;%
a            = 0;
e            = 0.0;
% U_speed = linspace(5,16,12);
U_speed = [0,linspace(3,15,13)];
h_sign  = -1; %-1 means EOM derived where airfoil going up in +h


%%

addpath('G:\My Drive\Codes\Aerodynamics-functions\Nonlinear\Bryant model\');

cprintf('*Keywords','[T multiplier]= %2.2f\n',T_multiplier)
cprintf('Strings','AR = %2d; PZT coverage ratio = %2.2f; PZT thickness ratio = %2.2f\n',AR,PZT_l_ratio,PZT_t_ratio)

aero_model    = 3;

% keyboard
%% Materials used
[Mat1] = read_material_data('Powerfilm_solar1.matl'); % Powerfilm_solar1.matl    QP16N_isotropic.matl
[Mat2] = read_material_data('QP16N_isotropic.matl'); %PTFE_ribbon_ANSYS.matl
Mats      = {Mat1, Mat2};


sim           = get_structural_simulation_options; %optional
% [seg, Mpoint] = geometry_layup_V2(Mats,L,ts,AR,PZT_thk_ratio,PZT_coverage,e,a);
% [seg, Mpoint] = geometry_layup_V3(Mats,L_s,thk_s,AR,PZT_thk,PZT_coverage,e,a);
[seg, Mpoint] = geometry_layup_V2(Mats,L_s,thk_s,AR,PZT_t_ratio,PZT_l_ratio,a,e);
%[seg, Mpoint]  = geometry_layup_V2(Mats,L_s,thk_s,AR,PZT_t_ratio,PZT_l_ratio,a,e); %Mats,L,ts,AR,PZT_thk_ratio,PZT_coverage,e)
% [seg, Mpoint] = geometry_layup_point_mass(Mats,L,ts,AR,PZT_thk_ratio,PZT_coverage,e);
segc          = segmental_discretizationV3(seg, Mats);     %(Mats,L,ts,AR,PZT_thk_ratio,PZT_coverage,e)
[panelc]      = panelwise_discretizationV3(segc,Mpoint);
Pcr           = pi^2*segc.EIbending(2)/(0.25*sum(segc.l)^2);

sim          = get_structural_simulation_options; %optional
Tactual      = Pcr*T_multiplier;
sim.T        = Tactual;
sim.nbending = 3;
sim.ntorsion = 2;
sim.theta0   = 0;
% keyboard
% sim.root_range_bending = [];

%% CURRENTLY WORKING ON THIS MODULE
% clc;
% [structural] = structural_module(segc,panelc,sim, Mpoint);%,matxRR  seg, panel, panelc, sim, Mpoint
FTM_size  = 1;
% idx_start = 76;
% idx_end      = idx_start - 1 + 235;
idx_start = 1;
idx_end      = numel(panelc.rhoA);
% [structural] = basic_structural_module(segc,panelc,sim, Mpoint,FTM_size ,idx_start,idx_end);  %,matxRR ,idx_start,idx_end
[structural,~,~,fHz_imbalance] = basic_structural_module(segc,panelc,sim, Mpoint,FTM_size ,idx_start,idx_end);

%% Electromechanical module
[theta_coupling0, Cp0]   = electromechanical_moduleV2(structural, panelc, segc); % seg, panel, panelc, sim, Mpoint , panelp ,panelp
structural.theta_coupL = reshape(theta_coupling0(1,1,:), [1,sim.nbending]) ;
structural.theta_coupR = reshape(theta_coupling0(1,3,:), [1,sim.nbending]);

structural.CpL         = Cp0(1);
structural.CpR         = Cp0(3);


%% Aeroelasticity module
Ropti                             = 22*1e3; %ohm,
[aero, sim]                       = get_aero_simulation_options(sim); %optional
sim.U_speed                       = U_speed;
[odedata,st,en,Uf_idx,Ropt]       = Bryant_aerodynamic_moduleV2(structural,panelc,sim,segc,aero,Ropti,h_sign);%,[1 1 1]
% keyboard
%% Post processing after aerolasticity solution
% h_sign = 1;
[tdomain, Udomain] = physical_output_postprocessing(odedata,st,en,sim,segc,panelc,structural,Ropt,h_sign);

%%
% clc
figt1    = time_domain_subplots(tdomain.t,tdomain.t_h,1000,0,0); %1000 = mm
figt2    = time_domain_subplots(tdomain.t,tdomain.t_alpha_eff,(180/pi),1,0);

%%
[Aswept,Power_aero] = aerodynamic_efficiency(odedata,st,en,structural.bending.PHI,structural.torsion.PSI,sim,panelc,Uf_idx,h_sign);
eta                 = Udomain.P(Uf_idx:end)./Power_aero; %Aswept,h_LE,h_TE,Power_aero
eta(isnan(eta))     = 0;

%%
for i = 1:numel(sim.U_speed)
    idx_limit   = find(odedata(i).ffx >= 250,1);
    fftx        = odedata(i).ffx(1:idx_limit);%fftx{i}
    ffty        = mag2db(odedata(i).ffy(1:idx_limit));%ffty{i}
    [~,idx_max] = max(ffty);
    fp(i)       = fftx(idx_max);
end
%% sending data out
% clc
data_out.A          = Udomain.A;
data_out.Ropt       = Ropt;
data_out.fp         = fp;
data_out.Uf_idx     = Uf_idx;
data_out.P_pzt      = Udomain.P;
data_out.theta      = rad2deg(Udomain.theta);    %deg
data_out.alpha_eff  = rad2deg(Udomain.alpha_eff); %deg
data_out.Aswept     = Aswept;
data_out.P_aero     = Power_aero;
data_out.eta        = eta.*100; %percentage
data_out.time       = tdomain.t;
data_out.t_theta    = tdomain.t_theta;
data_out.t_VL       = tdomain.t_VL;
data_out.Uspeed     = sim.U_speed;
data_out.PHI        = structural.bending.PHI;
data_out.PSI        = structural.torsion.PSI;
data_out.CthetaL    = structural.theta_coupL;
data_out.CthetaR    = structural.theta_coupR;
data_out.fheave     = structural.f_heave;
data_out.fpitch     = structural.f_pitch;
data_out.rmid       = panelc.rmid;
data_out.fHz_imbal  = fHz_imbalance;
data_out.Pcr        = Pcr;

% keyboard

%% Plotting for visualization
%%

figure
% subplot(6,1,1)
plot(sim.U_speed, data_out.A.*1000,'-b','LineWidth',2)
xlabel('U(m/s)')
ylabel('A(mm)')

figure
% subplot(6,1,2)
plot(sim.U_speed, data_out.P_pzt.*1000,'-b','LineWidth',2)
xlabel('U(m/s)')
ylabel('P(mW)')

figure
% subplot(6,1,3)
plot(sim.U_speed, data_out.alpha_eff,'-b','LineWidth',2)
hold on
plot(sim.U_speed, data_out.theta,'--r','LineWidth',2)
xlabel('U(m/s)')
ylabel('Angle(deg)')

figure
% subplot(6,1,4)
plot(sim.U_speed, data_out.eta,'-b','LineWidth',2)
xlabel('U(m/s)')
ylabel('\eta%')

figure
% subplot(6,1,5)
plot(sim.U_speed, data_out.fp,'-b','LineWidth',2)
xlabel('U(m/s)')
ylabel('f_p(Hz)')

figure
% subplot(6,1,6)
plot(sim.U_speed, data_out.Aswept.*1e6,'-b','LineWidth',2)
xlabel('U(m/s)')
ylabel('Aswept(mm^2)')

keyboard
% end