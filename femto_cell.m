clc; close all; clear all;

R = 1000;           % macro cell radius
rf = 20;            % femto cell radius
R0 = 1e3;           % out of cell interferance
SNR_targetdB = 10;  % SNR target in dB
rel = 0.95;         % reliability for SNR target
gam_mac = 3;        % pathloss exponent for the macro
gam_fem = 3;        % pathloss exponent for the femto
sig_sfdB = 8;       % shadowing std. dev (in dB)
p_ind = 0.5;        % prob femto user is indoorsa
p_act = 0.25;       % prob femto user is active
WdB = 10;           % wall loss
SNR_maxdB = 20;     % max avg SNR

N_exp = 10000;          % number of experiments
phi_f_values = [0 5 100 300]; % femto cell density values

F = 5;
kB = 1.38e-23;
T = 293;
B = 5e6;
NodBW = F+10*log10(kB*T*B);       % noise power

SINR_m = zeros(N_exp, length(phi_f_values)); % SINR for macro user 
C_m = zeros(N_exp, length(phi_f_values)); % macro user capacity
C_f_sum = zeros(N_exp, length(phi_f_values)); % femto user sum capacity

for p = 1:length(phi_f_values)
    phi_f = phi_f_values(p);
    PtMdB = getPtdB(R, NodBW, sig_sfdB, gam_mac, 0, WdB, rel, SNR_targetdB, 1e5);
    PtFdB = getPtdB(rf, NodBW, sig_sfdB, gam_fem, p_ind, WdB, rel, SNR_targetdB, 1e5); % transmitted power for femto BS      


    for i = 1:N_exp
        % get num. femto cells
        lam_f = phi_f*pi*(R/1000)^2;
        N_fem = max(1,(binornd(poissrnd(lam_f),p_act)));
            
        IfMdB = zeros(N_fem, 1); % interference from femto BS to macro user
        IfFdB = zeros(N_fem, 1); % interference from femto BS to femto user 
        SINR_f = zeros(1,N_fem); % SINR for femto user k
       
        % make macro cell
        macro_cell = struct();
        macro_cell.center_x = 0;
        macro_cell.center_y = 0;

        % drop one macro user randomly
        macro_cell.user_x = R * rand(1);
        macro_cell.user_y = R * rand(1);

        % drop phi_f femto cells somewhere in a radius R0=2e3
        femto_cells = struct();
        femto_users = struct();

        % create femto BS and users 
        for j = 1:N_fem
            % drop a femto cell
            femto_cell = randsphere(1, 2, R0);
            femto_cells(j).center_x = femto_cell(1);
            femto_cells(j).center_y = femto_cell(2);
        
            % drop one femto user for each cell
            femto_user = randsphere(1, 2, rf);
            femto_users(j).x = femto_user(1);
            femto_users(j).y = femto_user(2); 
        end

        % calculations for each femto BS and user
        for k = 1:N_fem
            % is femto user k inside or outside
            is_inside = binornd(1,1-p_ind,1);

            % distance calculations 
            d_fm = get_d(femto_cells(k).center_x, femto_cells(k).center_y, macro_cell.user_x, macro_cell.user_y); % between femto BS k and macro user
            d_ff = get_d(0,0, femto_users(k).x, femto_users(k).y); % between femto BS k and femto user k
            d_mf = get_d(0,0, femto_cells(k).center_x + femto_users(k).x, femto_cells(k).center_y + femto_users(k).y); % between macro BS and femto user k

            % power calculations
            PrFdB = get_PrdB(PtFdB, get_X(sig_sfdB), gam_fem, d_ff, get_H(), WdB*(1-is_inside), 1000); % recieved power for femto user k
            
            % interferance calculations
            IfMdB(k) = get_IdB(PtFdB, gam_mac, d_fm, get_X(sig_sfdB), WdB, get_H());% interference on macro user from femto cell k    
            ImFdB = get_IdB(PtMdB, gam_mac, d_mf, get_X(sig_sfdB), WdB*(1-is_inside), get_H()); % interferance on femto user k from macro BS 

            % calculate interferance from other femto BS to femto user k
            for l = 1:N_fem
                if(k ~= l)
                    d_ffL = get_d(femto_cells(l).center_x, femto_cells(l).center_y, femto_cells(k).center_x + femto_users(k).x, femto_cells(k).center_y + femto_users(k).y); % distance between femto BS l and femto user k
                    IfFdB(k) = get_IdB(PtFdB, gam_mac, d_fm, get_X(sig_sfdB), WdB*(1+is_inside), get_H()); % interference on macro user from femto cell k   
                end
            end
            SINR_f(k) = db2pow(PrFdB) / (db2pow(NodBW) + sum(db2pow(IfFdB))+sum( db2pow(ImFdB) ) );      
        end

        % distance between macro cell and macro user
        d_m = get_d(0, 0, macro_cell.user_x, macro_cell.user_y);

        % received power for macro user
        PrMdB = get_PrdB(PtMdB, get_X(sig_sfdB), gam_mac, d_m, get_H(), 0, SNR_maxdB+NodBW);

        SINR_m(i, p) = db2pow(PrMdB) / (db2pow(NodBW) + sum(db2pow(IfMdB)));
        
        % sum total femto user capacity
        count = 0;
        for m = 1:N_fem
             count = (log2(1 + SINR_f(m)) ) + count;
        end
        C_f_sum(i, p) = count;
    end

    % macro user capacity
    C_m(:, p) = (log2(1 + SINR_m(:, p)));
   
end

%Plotting
legend_labels = cell(length(phi_f_values), 1);
for p = 1:length(phi_f_values)
    phi_f = phi_f_values(p);
    legend_labels{2*p-1} = sprintf('phi_f = %d (C_m)', phi_f); 
    cdfplot(C_m(:, p)); hold on;
    if(phi_f ~= 0)
        legend_labels{2*p} = sprintf('phi_f = %d (C_f)', phi_f); 
        cdfplot(C_f_sum(:, p));
    end+
end

xlabel('Macro User Capacity (bps/Hz)');
ylabel('CDF');
legend(legend_labels, 'Location', 'Best');
xlim([0 10]);

% get transmitted power
function PtdB = getPtdB(r, NodBW, sig_sfdB, gam, p_ind, WdB, rel, SNR_targetdB, N_exp)

    %d = sqrt(rand(N_exp,1)*r^2);
    if (p_ind == 0)
        m_users = r * sqrt( rand(N_exp,1)) .* exp(1j*2*pi*rand(N_exp,1));
        PLdB = -10*gam*log10(abs(m_users))+sig_sfdB*randn(N_exp,1);
    else
        f_users = r * sqrt(rand(N_exp,1));
        PLdB = -10*gam*log10(abs(f_users))+sig_sfdB*randn(N_exp,1)-WdB*binornd(1, 1-p_ind, N_exp, 1);   
    end
    bins = -1000:0.001:1000;
    PLdB_cdf = cumsum(histc(PLdB,bins))/N_exp;
    [mincdf,n] = min(abs(PLdB_cdf-(1-rel)));
    PtdB = SNR_targetdB+NodBW-bins(n);
end

% euclidean distance between two points
function d = get_d(x1,y1,x2,y2)
    d = sqrt((x2 - x1)^2 + (y2 - y1)^2);
end

% get random point in sphere
function X = randsphere(m,n,r)
    X = randn(m, n);
    s2 = sum(X.^2, 2);
    X = X .* repmat(r * (gammainc(s2/2, n/2).^(1/n)) ./ sqrt(s2), 1, n);
end

% get the gaussian shadowing in dB
function XdB = get_X(sig_sfdB)
    XdB = sig_sfdB * randn();
end

% get Rayleigh fading. where abs(h^2)=H. 
function H = get_H()
    h = (randn() + 1i*randn()) / sqrt(2);
    H = abs(h)^2;
end

% get interferance
function IdB = get_IdB(PtdB, gam, d, XdB, WdB, H)
    IdB = PtdB+(-gam)*10*log10(d)+XdB-WdB+10*log10(H);
end

% get recieved power
function PrdB = get_PrdB(PtdB, XdB, gam, d, H, WdB, avgPr_max)
    PrdB = min(PtdB+XdB+(-gam)*10*log10(d)-WdB,avgPr_max)+pow2db(H);
end

