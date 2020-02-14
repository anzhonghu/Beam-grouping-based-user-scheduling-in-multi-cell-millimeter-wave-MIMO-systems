flag_kl1  = 0;
if 1==lll1
else
    for kc_c = 1 : Nr
        for ll_c = 1 : L
            if ll_c==l1
            else
                if K_ind(kc_c, ll_c) > 0
                    [~,ps] = max(beta_store(K_ind(kc_c, ll_c), (ll_c-1)*P+1:ll_c*P));
                    nx = floor((cos(thetab_store(K_ind(kc_c, ll_c), (ll_c-1)*P+ps)) * sin(phib_store(K_ind(kc_c, ll_c), (ll_c-1)*P+ps))+1)*naz*0.5);
                    ny = floor((sin(thetab_store(K_ind(kc_c, ll_c), (ll_c-1)*P+ps))+1)*nel*0.5);
                    nn = ny;
                    nnn=nx;
                    if block_store(k, (ll_c-1)*L+l1)>0
                    else
                        pp_c=1;
                        nx = floor((cos(thetabb_store(k, (ll_c-1)*L*P+(l1-1)*P+pp_c)) * sin(phibb_store(k, (ll_c-1)*L*P+(l1-1)*P+pp_c))+1)*naz*0.5);
                        ny = floor((sin(thetabb_store(k, (ll_c-1)*L*P+(l1-1)*P+pp_c))+1)*nel*0.5);
                        nn_c = ny;
                        nnn_c = nx;
                        if abs(nn_c-nn)+abs(nnn_c-nnn)<4
                            flag_kl1 = 1;
                            break;
                        else
                        end
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [~,ps] = max(beta_store(k, (l1-1)*P+1:l1*P));
                    nx = floor((cos(thetab_store(k, (l1-1)*P+ps)) * sin(phib_store(k, (l1-1)*P+ps))+1)*naz*0.5);
                    ny = floor((sin(thetab_store(k, (l1-1)*P+ps))+1)*nel*0.5);
                    nn = ny;
                    nnn=nx;
                    if block_store(K_ind(kc_c, ll_c), (l1-1)*L+ll_c)>0
                    else
                        pp_c=1;
                        nx = floor((cos(thetabb_store(K_ind(kc_c, ll_c), (l1-1)*L*P+(ll_c-1)*P+pp_c)) * sin(phibb_store(K_ind(kc_c, ll_c), (l1-1)*L*P+(ll_c-1)*P+pp_c))+1)*naz*0.5);
                        ny = floor((sin(thetabb_store(K_ind(kc_c, ll_c), (l1-1)*L*P+(ll_c-1)*P+pp_c))+1)*nel*0.5);
                        nn_c = ny;
                        nnn_c = nx;
                        if abs(nn_c-nn)+abs(nnn_c-nnn)<4
                            flag_kl1 = 1;
                            break;
                        else
                        end
                    end
                else
                end
            end
        end
    end
end