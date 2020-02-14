%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%%%capacity calculation%%%%%%%%%%%%%%%
capacity_temp_in = zeros(K,L);
for l1 = 1 : L%j, user cell
    nre=sum(C_sel(:, l1));
    if 0 == nre
        continue;
    else
    end
    for kc = 1 : nre
        k = K_ind(kc, l1);%%%%k
        if 1==kc
            [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
            nx = floor((cos(thetam_store(k, (l1-1)*P+ps)) * sin(phim_store(k, (l1-1)*P+ps))+1)*maz*0.5);
            ny = floor((sin(thetam_store(k, (l1-1)*P+ps))+1)*mel*0.5);
            nn = ny  * maz + nx + 1;
            wj = Um(:,nn);
        else
            [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
            nx = floor((cos(thetam_store(k, (l1-1)*P+ps)) * sin(phim_store(k, (l1-1)*P+ps))+1)*maz*0.5);
            ny = floor((sin(thetam_store(k, (l1-1)*P+ps))+1)*mel*0.5);
            nn = ny  * maz + nx + 1;
            wj_temp = [wj,Um(:,nn)];
            wj = wj_temp;
        end
    end
    interfandnoise = zeros(nre, 1);
    for l2 = 1 : L%%%l, bs cell
        nrel2=sum(C_sel(:, l2));
        if l2==l1 || 0== nrel2
            continue;
        else
            for kc = 1 : nrel2
                k = K_ind(kc, l2);
                if 1==kc
                    [~,ps] = max(beta_store(k, (l2-1)*P+1:l2*P));
                    nx = floor((cos(thetab_store(k, (l2-1)*P+ps)) * sin(phib_store(k, (l2-1)*P+ps))+1)*naz*0.5);
                    ny = floor((sin(thetab_store(k, (l2-1)*P+ps))+1)*nel*0.5);
                    nn = ny * naz + nx + 1;
                    Ul = U(:, nn);
                else
                    [~,ps] = max(beta_store(k, (l2-1)*P+1:l2*P));
                    nx = floor((cos(thetab_store(k, (l2-1)*P+ps)) * sin(phib_store(k, (l2-1)*P+ps))+1)*naz*0.5);
                    ny = floor((sin(thetab_store(k, (l2-1)*P+ps))+1)*nel*0.5);
                    nn = ny * naz + nx + 1;
                    Ul_temp = [Ul,U(:, nn)];
                    Ul = Ul_temp;
                end
            end
            Zjl = zeros(nre, nrel2);
            for kkk = 1 : nre
                k = K_ind(kkk, l1);%%%%k
                Zjl(kkk, :) = wj(:,kkk)' *  H(:,(l2-1)*K*L*n_arr+(l1-1)*K*n_arr+(k-1)*n_arr+1:(l2-1)*K*L*n_arr+(l1-1)*K*n_arr+k*n_arr) * Ul;
            end
            switch l2
                case 1
                    Zll = Zjj_1;
                case 2
                    Zll = Zjj_2;
                case 3
                    Zll = Zjj_3;
                otherwise
            end
            X = inv(Zll) * (inv(Zll))';
            for kc = 1 : nre
                interfandnoise(kc, 1) = interfandnoise(kc, 1) + (norm(Zjl(kc, :)/ Zll))^2 / trace(X);
            end
        end
    end
    interfandnoise = interfandnoise + 1/Pvsigma2;
    switch l1
        case 1
            Zjj = Zjj_1;
        case 2
            Zjj = Zjj_2;
        case 3
            Zjj = Zjj_3;
        otherwise
    end
    anotherf = trace(inv(Zjj) * (inv(Zjj))');
    for kc = 1 : nre
        k = K_ind(kc, l1);
        capacity_temp_in(k, l1) = log2(1+1 / abs(anotherf * interfandnoise(kc, 1)));
    end
end