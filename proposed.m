meth_ind = 4;
C_sel = zeros(K, L);
K_ind = zeros(Nr, L);
flagl1 = zeros(L, 1);
for lll1 = 1 : L
    l1=lll1;
    flagl1(l1, 1) = 1;
    hs = zeros(K, 1);
    gs = zeros(K, 1);
    interfusern = 0;
    for k = 1 : K
        cal_kl11;
        [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
        nx = floor((cos(thetam_store(k, (l1-1)*P+ps)) * sin(phim_store(k, (l1-1)*P+ps))+1)*maz*0.5);
        ny = floor((sin(thetam_store(k, (l1-1)*P+ps))+1)*mel*0.5);
        nn = ny  * maz + nx + 1;
        wj = Um(:,nn);
        nx = floor((cos(thetab_store(k, (l1-1)*P+ps)) * sin(phib_store(k, (l1-1)*P+ps))+1)*naz*0.5);
        ny = floor((sin(thetab_store(k, (l1-1)*P+ps))+1)*nel*0.5);
        nn = ny * naz + nx + 1;
        Ul = U(:, nn);
        if 1==flag_kl1
            interfusern = interfusern + 1;
            gs(k, 1) = abs(wj' * ...
                H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(k-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+k*n_arr) * Ul)^2;
        else
            hs(k, 1) = abs(wj' * ...
                H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(k-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+k*n_arr) * Ul)^2;
        end
    end
    [~, hs_ind] = sort(hs,'descend');
    [~, gs_ind] = sort(gs,'descend');
    hs_ind = [hs_ind(1:K-interfusern,1);gs_ind(1:interfusern,1)];
    kc = 0;
    for n = 1 : K
        k = hs_ind(n, 1);
        if 0==kc
            [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
            nx = floor((cos(thetam_store(k, (l1-1)*P+ps)) * sin(phim_store(k, (l1-1)*P+ps))+1)*maz*0.5);
            ny = floor((sin(thetam_store(k, (l1-1)*P+ps))+1)*mel*0.5);
            nn = ny  * maz + nx + 1;
            wj = Um(:,nn);
            nx = floor((cos(thetab_store(k, (l1-1)*P+ps)) * sin(phib_store(k, (l1-1)*P+ps))+1)*naz*0.5);
            ny = floor((sin(thetab_store(k, (l1-1)*P+ps))+1)*nel*0.5);
            nn = ny * naz + nx + 1;
            Ul = U(:, nn);
            kc = kc + 1;
            C_sel(k,l1) = 1;
            K_ind(kc, l1) = k;
        else
            [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
            nx = floor((cos(thetam_store(k, (l1-1)*P+ps)) * sin(phim_store(k, (l1-1)*P+ps))+1)*maz*0.5);
            ny = floor((sin(thetam_store(k, (l1-1)*P+ps))+1)*mel*0.5);
            nn = ny  * maz + nx + 1;
            wj_temp = [wj,Um(:,nn)];
            nx = floor((cos(thetab_store(k, (l1-1)*P+ps)) * sin(phib_store(k, (l1-1)*P+ps))+1)*naz*0.5);
            ny = floor((sin(thetab_store(k, (l1-1)*P+ps))+1)*nel*0.5);
            nn = ny * naz + nx + 1;
            Ul_temp = [Ul,U(:, nn)];
            Zll = zeros(kc+1,kc+1);
            for kkk = 1 : kc
                Zll(kkk, :) = wj_temp(:,kkk)' * ...
                    H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(K_ind(kkk, l1)-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+K_ind(kkk, l1)*n_arr)* Ul_temp;
            end
            kkk = kc + 1;
            Zll(kkk, :) = wj_temp(:,kkk)' * ...
                H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(k-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+k*n_arr) * Ul_temp;
            if cond(Zll)>1e3
            else
                wj = wj_temp;
                Ul = Ul_temp;
                kc = kc + 1;
                C_sel(k,l1) = 1;
                K_ind(kc, l1) = k;
            end
        end
        if kc == Nr
            break;
        else
        end
    end
    switch l1
        case 1
            kc = sum(C_sel(:,l1));
            Zjj_1 = zeros(kc,kc);
            for kkk = 1 : kc
                Zjj_1(kkk, :) = wj(:,kkk)' * ...
                    H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(K_ind(kkk, l1)-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+K_ind(kkk, l1)*n_arr)* Ul;
            end
        case 2
            kc = sum(C_sel(:,l1));
            Zjj_2 = zeros(kc,kc);
            for kkk = 1 : kc
                Zjj_2(kkk, :) = wj(:,kkk)' * ...
                    H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(K_ind(kkk, l1)-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+K_ind(kkk, l1)*n_arr)* Ul;
            end
        case 3
            kc = sum(C_sel(:,l1));
            Zjj_3 = zeros(kc,kc);
            for kkk = 1 : kc
                Zjj_3(kkk, :) = wj(:,kkk)' * ...
                    H(:,(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+(K_ind(kkk, l1)-1)*n_arr+1:(l1-1)*K*L*n_arr+(l1-1)*K*n_arr+K_ind(kkk, l1)*n_arr)* Ul;
            end
        otherwise
    end
end
capa_cal;
%%%%%%%%%%%%%%%%%%%%%%%%
capacity(variable_n, meth_ind) = capacity(variable_n, meth_ind) + sum(sum(capacity_temp_in)) * time_fre_resources;
jain_temp = 0;
for l1 = 1 : L
    for k = 1 : K
        jain_temp = jain_temp + capacity_temp_in(k, l1)^2;
    end
end
jain_temp = (sum(sum(capacity_temp_in)))^2 / K / L / jain_temp;
jain(variable_n, meth_ind) = jain(variable_n, meth_ind) + jain_temp;



