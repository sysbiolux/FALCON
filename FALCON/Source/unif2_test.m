function [U]=unif2_test(X)
%%% [U]=unif(X) computes the uniformity of the M by P matrix and returns
%%% the 1 by P vector of uniformities for each vector of M values.
%%% Uniformity is defined as the inverse of the average absolute deviation
%%% from the local density of values. Vectors of exactly equidistant values
%%% have infinite uniformity.
%%% This function considers all column vectors to be independent. To
%%% compute the uniformity of a P-dimensional matrix, the (weighted)
%%% average of the uniformities needs to be computed.
%%% Reference: De Landtsheer et al., Frontiers in Physiology, ...
%%% sebastien.delandtsheer@uni.lu
%%% sebdelandtsheer@gmail.com
%%% thomas.sauter@uni.lu


[T,S] = size(X);
U = 0;
for c = 1:(T-1)
    for cc = (c+1):T
        p1 = X(c, :);
        p2 = X(cc, :);
        Area = prod(abs(p1-p2));
        Inside_Idx = (X(:,2) > min(p1(2), p2(2))) & (X(:,1) > min(p1(1), p2(1))) & (X(:,1) < max(p1(1), p2(1))) & (X(:,2) < max(p1(2), p2(2)));
        Qtty = size(X(Inside_Idx,:),1);
        Expected = T;
        Density = (Qtty/Area);
        U = U + abs(Density - Expected);
    end
end

U = T ./U;





end