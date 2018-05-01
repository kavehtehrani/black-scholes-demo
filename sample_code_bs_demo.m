% Kaveh Tehrani
% Show casing greek surfaces for vanilla calls and puts

clc; close all; clear; tic;

%% INPUTS 
sigma = .2;
s0 = linspace(100*(1-sigma), 100*(1+sigma), 50);        % 1 std each way
ttm = linspace(1/252, 3/12, 50);                        % 1-day to 3-month option
k = 100;                                                % start at-the-money
rd = .03;
rf = 0.01;

%% CALCULATE THE OPTION VALUE + GREEKS
opt_prems = cell(length(s0), length(ttm));
greeks = cell(length(s0), length(ttm));
cp = struct();
for idx_s = 1:length(s0)
    for idx_t = 1:length(ttm)
        [ cp.c, cp.p ] = bs_option.calc_premium(s0(idx_s), k, rd, ttm(idx_t), sigma, rf);
        opt_prems{idx_s, idx_t} = cp;
        greeks{idx_s, idx_t} = bs_option.calc_greeks(s0(idx_s), k, rd, ttm(idx_t), sigma, rf);
    end
end

%% PLOT GREEK SURFACES
get_prem  = @(cp) cellfun(@(x) x.(cp), opt_prems);
get_greek = @(col, greek) cellfun(@(x) x.(greek)(col), greeks);

str_greeks = fields(greeks{1, 1});
[ xx, yy ] = meshgrid(ttm, s0');
for cp = 1:2
    % premium grid
    if cp == 1
        cur_prem = get_prem('c');
        col = 1;
        str_title = 'Call';
    else
        cur_prem = get_prem('p');
        col = 2;
        str_title = 'Put';
    end
    
    % greeks
    for idx_g = 1:length(str_greeks)
        cur_greek = get_greek(col, str_greeks{idx_g});
        
        figure;
        surf(xx, yy, cur_greek);
        xlabel('ttm'); ylabel('spot'); zlabel(strrep(str_greeks{idx_g}, '_', '\_'));
        title(upper([str_title ' ' strrep(str_greeks{idx_g}, '_', '\_')]))
    end
end

toc;
