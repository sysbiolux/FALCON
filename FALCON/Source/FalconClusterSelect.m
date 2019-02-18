function [ID] = FalconClusterSelect(params)
%%% FalconClusterSelect finds the best clustering for a set of values.
%%% 'params' is a vector, and ID is a vector of the indices of the clusters
%%% the algorithm performs 100 k-means clusterings with different values of
%%% k and reports for each one the one with the best silhouette. Then the
%%% optimal clustering is determined as the first maximum in the silhouette
%%% plot over all k.

%%% DEPRECATED as of V1.2 (February 2019)


N=length(params);
silA=zeros(1,N);
PC=zeros(N);
for k=1:N
    for k2=1:100
        idx=kmeans(params,k);
        score=mean(silhouette(params, idx));
        if score>silA(k), silA(k)=score; PC(k,:)=idx; end
    end
end
silA=[silA 0];
Pos_diff= find((silA(1:end-1)-silA(2:end))>0,1);

if ~isempty(Pos_diff)
    ID=(PC(Pos_diff,:))';
else
    ID=ones(length(params),1);
end

end