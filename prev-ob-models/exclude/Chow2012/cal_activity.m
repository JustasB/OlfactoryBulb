function [MC GC] = cal_activity(non_lin,CS,Wmg,Wgm,S,P,rm,rg)

    TEXT_OUT = 0;
    P_init = P;

    if non_lin == 0
        A = eye(size(S,1))+CS*Wmg*Wgm;
        MC = A\S;
        GC = Wgm*MC;
    end

    if non_lin == 1
        smooth = 0;
        P = zeros(size(S));
        stop_time = 20;
        for i = 1:size(S,2)
            tS = S(:,i);
            lastP = P_init(:,i);
            option = odeset('RelTol',1e-3,'AbsTol',1e-6,'OutputFcn',@out_fun2);
            [T,tempP] = ode23(@RHS2,[0 stop_time],P_init(:,i),option);
            P(:,i) = tempP(end,:);
            if TEXT_OUT == 1
                if T(end) >= stop_time
                    warning('max time exceeded!');
                end
            end
        end
        MC = P;
        GC = Wgm*rec(MC,rm,smooth);
    end
    
    function stop_ = out_fun2(t_,P_,flag_)
        stop_ = false;
        if strcmp(flag_,'init')
        elseif strcmp(flag_,'done')
        else
            if size(P_,2) > 1
                cond_ = norm(lastP-P_(:,end));
            else
                cond_ = norm(lastP-P_);
            end
            if cond_<1e-2
                stop_ = true;
            end
            if size(P_,2) > 1
                lastP = P_(:,end);
            else
                lastP = P_;
            end
        end  
    end % out_fun
    
    function dP = RHS2(t_,P_)
        dP = -P_ + tS - CS*Wmg*rec(Wgm*rec(P_,rm,smooth),rg,smooth);
    end % RHS

end % cal_activity