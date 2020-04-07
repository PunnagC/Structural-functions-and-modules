function [theta_seg, Cp_seg, theta, Cp] = electromechanical_module(structural,panel,segc,panelp)
% keyboard
IN_min = 3;
IN_max = 4;
narginchk(IN_min,IN_max);
nsegs  = numel(segc.l);
if nargin < IN_max
    [~,nbending] = size(structural.bending.PHI);
    theta        = zeros(2,nsegs,nbending);
    Cp           = zeros(2,nsegs);
    theta_seg    = theta;
    Cp_seg       = Cp;
else
    %     keyboard
    PHI_xx   = structural.bending.d2PHI_dx2;
    width    = panelp.chord;
    eP_T     = panelp.eP_T;
    eS       = panelp.eS;
    
    cumpnl   = [0, cumsum(segc.Npanel)]; %cumulative segmental panel numbers
    
    [~, cols_yu] = size(panelp.yu);%number of layers of piezos
    for i = 1:cols_yu
        for s = 1:nsegs
            st_idx       = cumpnl(s) + 1;
            en_idx       = cumpnl(s + 1);
            
            y_upper      = panelp.yu  (st_idx : en_idx, i); %checked ok
            y_lower      = panelp.yl  (st_idx : en_idx, i);%checked ok
            Eshi         = panelp.Eshi(st_idx : en_idx, i);%checked ok
            eP_Ti        = eP_T       (st_idx : en_idx, i); %checked ok
            widthi       = width      (st_idx : en_idx, i);
            eSi          = eS         (st_idx : en_idx, i);
            dr           = panel.dr   (st_idx : en_idx   );
            PHI_xx_s     = PHI_xx     (st_idx : en_idx, :);
            
            theta(i,s,:) = -(eP_Ti.*widthi.*Eshi.*dr.*(0.5*(y_upper.^2 - y_lower.^2)))'*PHI_xx_s;
            Cp   (i,s,:) = sum(eSi.*widthi.*Eshi.^2.*(y_upper - y_lower).*dr); %F, capacitance for the entire layer
            %             keyboard
        end
    end
    if cols_yu == 1 %unimorph
        theta_seg = theta;
        Cp_seg    = Cp;
    elseif cols_yu == 2 %unimorph
        theta_seg = sum(theta);
        Cp_seg    = sum(Cp);
    end
    
end


% keyboard
end