%% Black-Scholes european options valuation class (Garman and Kohlhagen specification)
% which is same as Black-Scholes for a continuously dividend paying stock)
% or alternatively european option on a foreign currency
% rd is domestic rate (or risk-free rate in BS)
% rf is foreign rate (or dividend rate (q) in BS)

classdef bs_option
    methods(Static = true)
        % d1 and d2 in BS
        function [ d1, d2 ] = calc_d(s0, k, rd, ttm, sigma, rf)
            bs_option.check_bs(s0, k, rd, ttm, sigma, rf);
            
            d1          = ( log(s0/k)+(rd-rf+0.5*(sigma^2))*ttm ) / (sigma*sqrt(ttm));
            d2          = d1 - sigma*sqrt(ttm);
        end
        
        % returns the BS option premium
        function [ call, put ] = calc_premium(s0, k, rd, ttm, sigma, rf)
            if ttm == 0
                % at expiry
                call        = max(0, s0-k);
                put         = max(0, k-s0);
            else
                % otherwise
                [d1, d2]    = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
                
                call        = s0*exp(-rf*ttm)*normcdf(d1) - k*exp(-rd*ttm)*normcdf(d2);
                put         = k*exp(-rd*ttm)*normcdf(-d2) - s0*exp(-rf*ttm)*normcdf(-d1);
            end
        end
        
        % returns the BS implied volatility that matches the given price
        function iv = calc_implied_vol(opt_type, prem, s0, k, rd, ttm, rf)
            init_iv_guess       = .2;
            max_iv              = 2;
            
            % solve within condition of 0 < iv < 2
            func_iv             = @(sigma)          nested_val(sigma, opt_type, prem, s0, k, rd, ttm, rf);
            [ iv, resnorm, res, exitflag ]          = lsqnonlin(func_iv, init_iv_guess, 0, max_iv);
            if exitflag ~= 1 && exitflag ~= 3,      iv = [];    end
            
            % nested function to value the call or put for implied volatility solution
            function val = nested_val(sigma, opt_type, prem, s0, k, rd, ttm, rf)
                [ c, p ]                            = bs_option.calc_premium(s0, k, rd, ttm, sigma, rf);
                if strcmpi(opt_type, 'call') || strcmpi(opt_type, 'c')
                    val = c - prem;
                elseif strcmpi(opt_type, 'put') || strcmpi(opt_type, 'p')
                    val = p - prem;
                else                                error('unknown option type encountered.');
                end
            end
        end
        
        % returns the strike that matches specified delta
        function k = calc_strike_by_delta(opt_type, tgt_delta, s0, rd, ttm, sigma, rf, b_display)
            UPPER_X                                 = 5;
            if abs(tgt_delta) > 1
                error('abs(delta) cannot be greater than one.');        
            end
            
            init_k_guess                            = s0;
            
            func_k              = @(k)              nested_delta(opt_type, tgt_delta, s0, k, rd, ttm, sigma, rf);

            
            [ k, resnorm, res, exitflag ]           = lsqnonlin(func_k, init_k_guess, 0, s0*UPPER_X);
            if exitflag ~= 1 && exitflag ~= 3,      
                k = [];     
            end

            function val = nested_delta(opt_type, tgt_delta, s0, k, rd, ttm, sigma, rf)
                delta                               = bs_option.calc_delta(s0, k, rd, ttm, sigma, rf);
                if strcmpi(opt_type, 'call') || strcmpi(opt_type, 'c')
                    val = tgt_delta - delta(1);
                elseif strcmpi(opt_type, 'put') || strcmpi(opt_type, 'p')
                    val = tgt_delta - delta(2);
                else                                error('unknown option type encountered.');
                end
            end
        end

        % option greeks
        function delta = calc_delta(s0, k, rd, ttm, sigma, rf)
            [d1, ~]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            delta       = [ exp(-rf*ttm)*normcdf(d1) exp(-rf*ttm)*(normcdf(d1)-1) ];
        end
        
        function gamma = calc_gamma(s0, k, rd, ttm, sigma, rf)
            [d1, ~]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            gamma       = [ (normpdf(d1)*exp(-rf*ttm))/(s0*sigma*sqrt(ttm)) (normpdf(d1)*exp(-rf*ttm))/(s0*sigma*sqrt(ttm)) ];
        end
        
        function theta = calc_theta(s0, k, rd, ttm, sigma, rf)
            [d1, d2]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            theta       = [ -(s0*exp(-rf*ttm)*normpdf(d1)*sigma)/(2*sqrt(ttm)) + rf*s0*exp(-rf*ttm)*normcdf(d1) - rd*k*exp(-rd*ttm)*normcdf(d2) ...
                -(s0*exp(-rf*ttm)*normpdf(d1)*sigma)/(2*sqrt(ttm)) - rf*s0*exp(-rf*ttm)*normcdf(-d1) + rd*k*exp(-rd*ttm)*normcdf(-d2) ];
        end
        
        function vega = calc_vega(s0, k, rd, ttm, sigma, rf)
            [d1, ~]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            vega        = [ s0*exp(-rf*ttm)*normpdf(d1)*sqrt(ttm) s0*exp(-rf*ttm)*normpdf(d1)*sqrt(ttm) ];
        end
        
        % also known as domestic rho for fx options
        function rho = calc_rho(s0, k, rd, ttm, sigma, rf)
            [~, d2]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            rho         = [ k*ttm*exp(-rd*ttm)*normcdf(d2) -k*ttm*exp(-rd*ttm)*normcdf(-d2) ];
        end
        
        % also known as foreign rho for fx options
        function q_rho = calc_q_rho(s0, k, rd, ttm, sigma, rf)
            [d1, ~]     = bs_option.calc_d(s0, k, rd, ttm, sigma, rf);
            q_rho       = [ -s0*ttm*exp(-rf*ttm)*normcdf(d1) s0*ttm*exp(-rf*ttm)*normcdf(-d1) ];
        end
        
        % packages all greeks in one function
        function [ greeks ] = calc_greeks(s0, k, rd, ttm, sigma, rf)
            greeks.delta    = bs_option.calc_delta(s0, k, rd, ttm, sigma, rf);
            greeks.gamma    = bs_option.calc_gamma(s0, k, rd, ttm, sigma, rf);
            greeks.theta    = bs_option.calc_theta(s0, k, rd, ttm, sigma, rf);
            greeks.vega     = bs_option.calc_vega(s0, k, rd, ttm, sigma, rf);
            greeks.rho      = bs_option.calc_rho(s0, k, rd, ttm, sigma, rf);
            greeks.q_rho    = bs_option.calc_q_rho(s0, k, rd, ttm, sigma, rf);
        end
        
        % fx option quotation convention, vol is quoted for delta-neutral spot straddles
        function k_atm = calc_atm_strike(s0, rd, ttm, sigma, rf)
            bs_option.check_bs(s0, 0, rd, ttm, sigma, rf);
            
            k_atm       = s0*exp((rd-rf+.5*(sigma^2))*ttm);
        end
        
        % translates quoted vol and delta to a strike
        function k = calc_strike(s0, rd, ttm, sigma, rf, delta, b_fwd_delta)
            if nargin < 7
                b_fwd_delta = true;
            end
            
            bs_option.check_bs(s0, 0, rd, ttm, sigma, rf);
            if abs(delta) > 1
                error('bsoption: delta absolute value cannot be more than one.');
            end
            
            if b_fwd_delta
                mult    = 1;
            else
                mult    = exp(rf*ttm);
            end
            
            k_call      = s0*exp((rd-rf)*ttm) * exp( -sigma*sqrt(ttm)*norminv(mult*delta) + .5*(sigma^2)*ttm  );
            k_put       = s0*exp((rd-rf)*ttm) * exp( sigma*sqrt(ttm)*norminv(-mult*delta) + .5*(sigma^2)*ttm  );
            
            if delta < 0,       k = k_put;
            else                k = k_call;
            end
        end
        
        % translate deltas into call deltas (absolute value forward put and call deltas sum to 1)
        function call_delta = harmonize_delta(delta)
            call_delta = delta;
            if any(abs(delta) > 1)
                error('delta absolute value cannot be more than 1');
            else
                neg_rng = delta < 0;
                call_delta(neg_rng)     = 1 + delta(neg_rng);
                % atm dns has call delta of .5 and put delta of -.5
                call_delta(delta == 0)  = .5;
            end
        end
    end
    
    methods (Access = protected, Static = true)
        % all paramters need to be non-negative except the rates
        function b_valid = check_bs(s0, k, rd, ttm, sigma, rf)
            if s0 < 0
                error('bscheck: negative spot price.');
            elseif k < 0
                error('bscheck: negative strike price.');
            elseif ttm < 0
                error('bscheck: negative time to maturity.');
            elseif sigma < 0
                error('bscheck: negative volatility.');
            elseif rf < 0
                warning('bs_option:bscheck', 'negative asset yield');
                b_valid = 1;
            elseif rd < 0
                warning('bs_option:bscheck', 'negative risk free rate');
                b_valid = 1;
            else
                b_valid = 1;
            end
        end
    end
end

