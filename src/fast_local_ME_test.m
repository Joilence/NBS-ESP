function [pv_left_mat,odds_mat] = ...
    fast_local_ME_test(seed_gene_set,p2m,network,ngh_type,cancer_type)

ngene = size(p2m,2);
pv_left_mat = ones(ngene,ngene);
odds_mat = zeros(ngene,ngene);
tct = 0;
s = sum(p2m);
for g= seed_gene_set
    tct = tct + 1;
    if mod(tct,1000)==0
        fprintf('%s finsihed %f\n',cancer_type,tct/length(seed_gene_set));
    end
    if ngh_type==1
        ngh = find(network(g,:)>0);
    else
        ngh = seed_gene_set;
    end
    gx = p2m(:,g)';
    for ng = ngh
        if s(g)<s(ng)
            continue;
        end
        px = p2m(:,ng)';
        if sum(px)==0 || sum(gx)==0
            continue
        end
        ct = crosstab(px,gx);
        if min(size(ct))<2
            continue;
        end
        
        [~,pv_left,left_stats] = fishertest(ct,'Tail','left');
        pv_left_mat(g,ng) = pv_left;
        pv_left_mat(ng,g) = pv_left;
        
        odds_mat(g,ng) = left_stats.OddsRatio;
        
        odds_mat(ng,g) = left_stats.OddsRatio;
    end
end


end

