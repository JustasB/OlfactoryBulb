function main
    neurogenesis(0);
end % main 

function objfun = neurogenesis(run_num)

% parameters
    % output
        TEXT_OUT = 1;
        PLOT_OUT = 1; ANIME_OUT = 1;
        PLOT2D_OUT = 1; ANIME2D_OUT = 0;
    % stimulus
        sim = 40;
        odor_names = 'limonene(+)_ster limonene(_)_ster propylpropionate_es3 ethylbutyrate_es3 isopropylbenzene_ModuleC1 cyclohexanone_SG18 acetone methylacetate_SG19 cycloheptanelow_cycloalk propanol_simp_2500 isoamylbutyrate_est1 butyricacid_aci1 hexanal_ald1 ethylbenzene_HC';
        choose = 1:14;
    % network
        non_lin = 0;
        conn = 4;
        CS = 0.002;
    % survival
        ts = 0.1;
        gamma = 5/ts;
        th = 0.2;
        rm = 0;
        rg = 0;
    % stepping
        cont_density = 0;
        exp_time = 5000; step = round(exp_time/100); dt = 5;
    % tracking vs T
        tracking = 1; % 0: none, 1: in S
        track_pairs = [1 2; 3 4];
        
    % ...
        S0 = 1; Sstr = 1; MC_per_Glom = 1;
        TYPE = 1; % 1: pearson 2: L2
        prob_conn = 0; % only for cont_density = 0
        Na = 1;
        perm_ratio = 0; % only for cont_density = 0
        minv = 0; maxv = 1;
        
% setup stimulus
        
% setup stimulus
    if TEXT_OUT == 1
        fprintf('\n--- Start ---\n   - initializing odor\n');
    end
    [Sall, coord, metric, name] = gara(14, MC_per_Glom, odor_names, sim);
    Nc = size(Sall,1);
    Ns = length(choose);
    S = S0+Sstr*(Sall(:,choose));
    S_name = cell(1,size(S,2));
    for i = 1:size(S,2)
        S_name{i} = name{choose(i)};
    end
    
    if TEXT_OUT == 1
        fprintf('       - Nc = %d\n', Nc);
        fprintf('           corr = %f\n', mean_excluNaN(uptri_1d(corr(S,TYPE))));
    end
    if PLOT_OUT == 1;
        setup_Pplot(S,corr(S,TYPE),corr(S',TYPE),rg,1);
        drawnow;
    end
    if PLOT2D_OUT == 1;
        setup_Pplot2D(S,coord,101,S_name);
        drawnow;
    end
    
% setup network
    if TEXT_OUT == 1
        fprintf('   - initializing network\n');
    end
    if cont_density == 1
        option = odeset('Stats','off','RelTol',1e-3,'AbsTol',1e-8);
        Isize = nchoosek(Nc,conn);
        perm = nchoosek(1:Nc,conn);
        C = zeros(Nc,Isize);
        for i = 1:Isize
            for c = 1:conn
                C(perm(i,c),i) = 1;
            end
        end
        N = rand(Isize,1);
        Wmg = C*diag(N); Wgm = C';
        Iage = ones(1,Isize);
        Imark = zeros(1,Isize);
    else
        Isize = 2*Nc; N = ones(Isize,1);
        Wmg = zeros(Nc, Isize); Wgm = zeros(Isize, Nc);
        Iage = -ones(1,Isize);
        Imark = zeros(1,Isize);  
    end
    [P, I] = cal_activity(non_lin,CS,Wmg,Wgm,S,S,rm,rg);
    time_axis = 0:step:exp_time;
    N_t = NaN*ones(length(time_axis),Isize);
    Pcorr_t = NaN*ones(1,length(time_axis));
    Tcorr_t = NaN*ones(1,length(time_axis));
    Pangle_t = NaN*ones(1,length(time_axis));
    Tangle_t = NaN*ones(1,length(time_axis));
    Pfoc_t = NaN*ones(1,length(time_axis));
    CV_t = NaN*ones(1,length(time_axis));
    CVid_t = NaN*ones(Ns,length(time_axis));
    F_t = NaN*ones(Ns,Ns,length(time_axis));
    if cont_density == 1
        N_t(1,:) = N;
    end
    if PLOT_OUT == 1;
        HP1d = setup_Pplot(P,corr(P,TYPE),corr(P',TYPE),rg,2,Wmg,Wgm,time_axis,N_t);
        [HI, Iaxis] = setup_Iplot(cont_density,time_axis,I,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(1,:),3);
        VST = cell(1,2);
        VST{1} = Pcorr_t; VST{2} = Tcorr_t;
        line_style = cell(1,2);
        line_style{1} = '-b'; line_style{2} = '-r';
        line_name = cell(1,2);
        line_name{1} = Pcorr_t; line_name{2} = Tcorr_t;
        Hinfo = setup_Infoplot(time_axis,VST,line_style,line_name,corr(S,TYPE),corr(P,TYPE),4);
        drawnow;
    end
    if PLOT2D_OUT == 1;
        setup_Pplot2D(P,coord,102,S_name);
        drawnow;
    end
    
% step
    if TEXT_OUT == 1
        fprintf('   - running\n');
    end
    Pcorr_t(1) = mean_excluNaN(uptri_1d(corr(P,TYPE)));
    Pangle_t(1) = mean_excluNaN(uptri_1d(corr_angle(corr(S,TYPE),corr(P,TYPE))));
    Pfoc_t(1) = mean_excluNaN(focality(P,metric));
    CV_t(1) = std(mean(P,2))/mean(mean(P,2));
    CVid_t(:,1) = std(P)./mean(P);
    F_t(:,:,1) = corr(P);
    if tracking == 1
        Tcorr_t(1) = mean_excluNaN(cal_track_corr(track_pairs,P));
        Tangle_t(1) = mean_excluNaN(corr_angle(cal_track_corr(track_pairs,S),cal_track_corr(track_pairs,P)));
    end
    if TEXT_OUT == 1
        fprintf('           corr = %f\n', Pcorr_t(1));
    end
    for i = 1:round(exp_time/step)
        if TEXT_OUT == 1
            fprintf('       - run num = %d, time = %f\n',run_num,i*step);
        end
        
        if cont_density == 1
            [ignore,N] = ode23(@RHS,[0 step],N_t(i,:),option);
            N_t(i+1,:) = N(end,:);
            Wmg = C*diag(N(end,:)); Wgm = C';
        else
            for j = 1:round(step/dt)
                add_cell;
                [P, I] = cal_activity(non_lin,CS,Wmg,Wgm,S,P,rm,rg);
                remove_cell;
            end
        end
        [P, I] = cal_activity(non_lin,CS,Wmg,Wgm,S,P,rm,rg);
        Pcorr_t(i+1) = mean_excluNaN(uptri_1d(corr(P,TYPE)));
        Pangle_t(i+1) = mean_excluNaN(uptri_1d(corr_angle(corr(S,TYPE),corr(P,TYPE))));
        Pfoc_t(i+1) = mean_excluNaN(focality(P,metric));
        CV_t(i+1) = std(mean(P,2))/mean(mean(P,2));
        CVid_t(:,i+1) = std(P)./mean(P);
        F_t(:,:,i+1) = corr(P);
        if tracking == 1
            Tcorr_t(i+1) = mean_excluNaN(cal_track_corr(track_pairs,P));
            Tangle_t(i+1) = mean_excluNaN(corr_angle(cal_track_corr(track_pairs,S),cal_track_corr(track_pairs,P)));
        end
        
        if TEXT_OUT == 1
            fprintf('           corr = %f\n', Pcorr_t(i+1));
        end

        if ANIME_OUT == 1;
            update_Pplot(P,corr(P,TYPE),corr(P',TYPE),rg,HP1d,Wmg,Wgm,N_t);
            update_Iplot(cont_density,time_axis,I,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(i,:),HI);
            VST = cell(1,2);
            VST{1} = Pcorr_t; VST{2} = Tcorr_t;
            update_Infoplot(VST,corr(P,TYPE),Hinfo);
            drawnow;
        end
        if ANIME2D_OUT == 1;
            setup_Pplot2D(P,coord,102,S_name);
            drawnow;
        end
    end
    
% end

    if PLOT_OUT == 1;
        update_Pplot(P,corr(P,TYPE),corr(P',TYPE),rg,HP1d,Wmg,Wgm,N_t);
        update_Iplot(cont_density,time_axis,I,corr(I(Iage>=0,:),TYPE),Iage,Wmg,N_t(i,:),HI);
        VST = cell(1,2);
        VST{1} = Pcorr_t; VST{2} = Tcorr_t;
        update_Infoplot(VST,corr(P,TYPE),Hinfo);
        drawnow;
    end
    if PLOT2D_OUT == 1;
        setup_Pplot2D(P,coord,102,S_name);
        drawnow;
    end
    
    objfun = return_val;
    if TEXT_OUT == 1
        fprintf('--- End ---\n');
    end

% nested function definition
    function dN_ = RHS(ignore,N1_)
        N_ = N1_;
        Wmg = C*diag(N_); Wgm = C';
        [ignore,G_] = cal_activity(non_lin,CS,Wmg,Wgm,S,P,rm,rg);
        P_ = survival(G_);
        P_ = P_ + 1e-10*randn(size(P_));
        B_ = Na*Nc/Isize/conn;
        dN_ = B_ + log(P_).*N_;
    end % RHS

    function prob_ = survival(G_)
        Ca_ = sum(rec(G_,th,0),2);
        prob_ = (tanh((Ca_-ts)*gamma)+1)*(maxv-minv)/2+minv;
    end % survival
    
    function val = remove_cell
        prob_ = survival(I);
        surv_ = floor(prob_ + rand(size(prob_)));
        IX_ = find((surv_==0)&(Iage'>=0));
        Iage(IX_) = -1;
        Imark(IX_) = 0;
        Wmg(:,IX_) = 0;
        Wgm(IX_,:) = 0;
        val = length(IX_);
    end % remove_cell
    
    function val = add_cell
        IX_ = find(Iage>=0);
        Iage(IX_) = Iage(IX_)+dt;
        IX_ = find(Iage<0);
        if length(IX_) < round(dt*Na*Nc)
            add_space(round(dt*Na*Nc)-length(IX_));
            IX_ = find(Iage<0);
        end
        if prob_conn == 0
            temp_ = [ones(conn,1); zeros(Nc-conn,1)];
            for i_ = 1:round(dt*Na*Nc)
                Iage(IX_(i_)) = 0;
                temp_ = temp_(randperm(Nc));
                Wmg(:, IX_(i_)) = temp_;
                IX2_ = randperm(length(temp_));
                temp1_ = temp_(IX2_(1:round(perm_ratio*length(temp_))));
                temp2_ = temp_(IX2_(1+round(perm_ratio*length(temp_)):end));
                Wgm(IX_(i_),IX2_) = [temp1_(randperm(length(temp1_))); temp2_]';
            end
        else
            conn_prob_ = conn/Nc;
            for i_ = 1:round(dt*Na*Nc)
                Iage(IX_(i_)) = 0;
                if (marking == 1) && (i*step >= marking_t(1)) && (i*step < marking_t(2))
                    Imark(IX_(i_)) = 1;
                end
                temp_ = floor(rand(Nc,1)+conn_prob_);
                Wmg(:, IX_(i_)) = temp_;
                IX2_ = randperm(length(temp_));
                temp1_ = temp_(IX2_(1:round(perm_ratio*length(temp_))));
                temp2_ = temp_(IX2_(1+round(perm_ratio*length(temp_)):end));
                Wgm(IX_(i_),IX2_) = [temp1_(randperm(length(temp1_))); temp2_]';
            end
        end
        val = round(dt*Na*Nc);
    end % add_cell
    
    function add_space(short)
        previous_size = Isize;
        tempI = I;
        tempIage = Iage;
        tempImark = Imark;
        tempWmg = Wmg;
        tempWgm = Wgm;
        
        Isize = Isize + 5*short;
        I = zeros(Isize, Ns);
        Iage = -ones(1, Isize);
        Imark = zeros(1, Isize);
        Wmg = zeros(Nc, Isize);
        Wgm = zeros(Isize, Nc);
        
        I(1:previous_size, :) = tempI;
        Iage(1, 1:previous_size) = tempIage;
        Imark(1, 1:previous_size) = tempImark;
        Wmg(:, 1:previous_size) = tempWmg;
        Wgm(1:previous_size,:) = tempWgm;
                
        if ANIME_OUT==1
            set(Iaxis, 'YLim', [1 Isize]);
        end
    end % add_space
    
    function val = cal_track_corr(track_pairs_,S_)
        corr_ = zeros(1,size(track_pairs_,1));
        for i_ = 1:size(track_pairs_,1)
            temp_ = corr(S_(:,track_pairs_(i_,:)),TYPE);
            corr_(i_) = temp_(1,2);
        end
        val = corr_;
    end % cal_track_corr
    
    function val = return_val
        
        val = 0;

    end % return_val

end % neurogenesis
