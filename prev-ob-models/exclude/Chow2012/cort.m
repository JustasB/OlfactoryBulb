function cort

    all_Ns = 33;
    MC_per_Glom = 1;
    sim = 20;
    odor_names = 'limonene(+)_ster limonene(_)_ster propylpropionate_es3 ethylbutyrate_es3 isopropylbenzene_ModuleC1 cyclohexanone_SG18 acetone methylacetate_SG19 cycloheptanelow_cycloalk propanol_simp_2500 isoamylbutyrate_est1 butyricacid_aci1 hexanal_ald1 ethylbenzene_HC D_carvone_ster L_carvone_ster 1_butanol_SG12 1_heptanol_OH2 1_hexanol_alc2 aceticacid_aci1 caproicacid_aci1_7.2 dodecanal_SG20 valericacid_aci1 2_hexanone_fgrp_250 amylacetate_es3 butylbutyrate_banana1 citronellol_C10Alcs eucalyptol_polycy me_salicylate_SG8 menthylisoval_mint_1_4 pentanol_fgrp_250 terp4ol(+)_ster terp4ol(_)_ster';
    [Sall, ~, ~, ~] = gara(all_Ns, MC_per_Glom, odor_names, sim);
    Sall = 1+Sall;
    Nc = size(Sall,1);
    
    steps = 250;
    Isize = 2500;
    conn = 8;
    CS = 0.02;
    rate = 0.01;
    
    choose = 1:8;
    Ns = length(choose);
    S = Sall(:,choose);
    
    Wmg = zeros(Nc, Isize);
    G_cort = zeros(1,Isize);
    nd = round(Isize*rate);
    
    [P, ~] = cal_activity(0,CS,Wmg,Wmg',S,S);
    P(P<0) = 0;
    [Pall, ~] = cal_activity(0,CS,Wmg,Wmg',Sall,Sall);
    Pall(Pall<0) = 0;
    for i = 1:steps
        CP = C(P);
        if i ~=1
            cort_act = CP;
            Gsum = cort_act(G_cort);
            rand_str = max(abs(Gsum));
            [~, IX] = sort(Gsum+rand_str*randn(size(Gsum)));
            IX = IX(1:nd);
        else
            IX = 1:Isize;
        end
        prob = [0 cumsum(CP/sum(CP))];
        prob_recp = [zeros(1,all_Ns); cumsum(Pall./repmat(sum(Pall),Nc,1))];
        for j = 1:length(IX)
            seed = rand; temp = find((prob-seed)<=0); temp = temp(end);
            G_cort(IX(j)) = temp;
            prob_temp = prob_recp(:,G_cort(IX(j)));
            k = 0;
            Wmg(:,IX(j)) = 0;
            while k < conn
                seed = rand; temp = find((prob_temp-seed)<=0); temp = temp(end);
                if Wmg(temp,IX(j))~=1
                    Wmg(temp,IX(j)) = 1;
                    k = k+1;
                end
            end
        end
        [P, ~] = cal_activity(0,CS,Wmg,Wmg',S,S);
        P(P<0) = 0;
        [Pall, ~] = cal_activity(0,CS,Wmg,Wmg',Sall,Sall);
        Pall(Pall<0) = 0;

        figure(1);
        subplot(2,2,1);
        imagesc(Wmg*Wmg'-diag(NaN*ones(1,Nc)));
        subplot(2,2,3);
        hist(G_cort,0:all_Ns);
        xlim([0.5 all_Ns+0.5]);
        subplot(2,2,2);
        imagesc(Pall);
        subplot(2,2,4);
        imagesc(corrcoef(Pall)); caxis([-1 1]);
        
        drawnow;
    end
    
    function val = C(p_)
        val = zeros(1,all_Ns);
        for j_ = 1:size(p_,2)
            for i_ = 1:all_Ns
                temp_ = corrcoef([Pall(:,i_) p_(:,j_)]);
                val(i_) = val(i_)+max(0,temp_(1,2));
            end
        end
    end

end