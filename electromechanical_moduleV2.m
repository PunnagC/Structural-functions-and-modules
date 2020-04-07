function [theta_seg, Cp_seg, theta, Cp] = electromechanical_moduleV2(structural,panelc,segc)
% keyboard
IN_min = 2;
IN_max = 3;
narginchk(IN_min,IN_max);
% ns    = numel(segc.l);
[nl, ns] = size(segc.layer_arrangement); %number of structural segments of the blade

if nargin < IN_max
    [~,nbending] = size(structural.bending.PHI);
    theta        = zeros(nl,ns,nbending);
    Cp           = zeros(nl,ns);
    theta_seg    = theta;
    Cp_seg       = Cp;
else
    
% keyboard
PHI_xx   = structural.bending.d2PHI_dx2;
width    = panelc.chord;
eP_T     = panelc.eP_T;
eS       = panelc.eS;

cumpnl   = [0, cumsum(segc.Npanel)]; %cumulative segmental panel numbers

for i = 1:ns     %sweeping through segments

    [nl]    = nnz(segc.layer_arrangement(:,i)); %number of structural segments of the blade
    st_idx  = cumpnl(i) + 1;
    en_idx  = cumpnl(i + 1);

    widthi   = width    (st_idx : en_idx);   %npx1
    dr       = panelc.dr(st_idx : en_idx);   %npx1
    PHI_xx_i = PHI_xx   (st_idx : en_idx, :);%npxnh
    %         disp(['Segment=',num2str(i)])
    
    for j = 1:nl %sweeping through layers j = 1 (bottom most layer)
        %             disp(['Layer=',num2str(j)])
        %                   keyboard
    
        y_upper      = panelc.yu{i}( :, j) ; %npx1 checked ok ( :, j)
        y_lower      = panelc.yl{i}( :, j) ; %npx1 %checked ok ( :, j)
        
        Eshi         = panelc.Eshi{i}( :, j); %npx1 %checked ok
        eP_Ti        = eP_T       {i}( :, j); %npx1 %checked ok
        eSi          = eS         {i}( :, j); %npx1
%         keyboard
        theta_modal  = -(eP_Ti.*widthi.*dr.*(0.5.*Eshi.*(y_upper.^2 - y_lower.^2)))'*PHI_xx_i; %z axis of 3D matrix is the mode number
        thetaji(j,i) =  sum(theta_modal);
        
        theta(j,i,:) = theta_modal;
        Cp   (j,i,:) = sum(eSi.*widthi.*Eshi.^2.*abs(y_upper - y_lower).*dr); %F, capacitance for the entire layer
%                      keyboard
    end

end
% keyboard
% theta = sum(thetaji);
%% Hardcoding 2nd mode coupling coefficients
[rx,cy,mz] = size(theta);
% keyboard
% if cy >= 3 & mz >= 2
%     theta(:,3,2) = -theta(:,3,2); % such that there is not cancellation
% end

% for i = 1:ns
%     for j = 1:nl
    %%
%     theta
%     Cp
% keyboard
    if j > 2 %bimorph
        theta_seg = sum(theta);
%         theta_seg = sum(thetaji);
        Cp_seg    = sum(Cp);
    else
        theta_seg = theta;
        Cp_seg    = Cp;
    end
    
    
    
    
%     if cols_yu == 1 %unimorph
%         theta_seg = theta;
%         Cp_seg    = Cp;
%     elseif j > 2 %bimorph
%         theta_seg = sum(theta);
%         Cp_seg    = sum(Cp);
%     end
    
end
% keyboard
end