capacity_cumu = zeros(K, L);
meth_ind = 2;
for tfr = 1 : time_fre_resources
    C_sel = zeros(K, L);
    K_ind = zeros(Nr, L);
    for l1 = 1 : L
        [~, Csel] = sort(capacity_cumu(:, l1), 'ascend');
        kc = 0;
        for n = 1 : K
            if 1==n
                k = Csel(n,1);
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
                k = Csel(n,1);
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
                if cond(Zll * Zll')>1e3
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
    capacity_cumu = capacity_cumu + capacity_temp_in;
end
capacity(variable_n, meth_ind) = capacity(variable_n, meth_ind) + sum(sum(capacity_cumu));
jain_temp = 0;
for l1 = 1 : L
    for k = 1 : K
        jain_temp = jain_temp + capacity_cumu(k, l1)^2;
    end
end
jain_temp = (sum(sum(capacity_cumu)))^2 / K / L / jain_temp;
jain(variable_n, meth_ind) = jain(variable_n, meth_ind) + jain_temp;



