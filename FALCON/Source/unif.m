function [U]=unif(X)
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
X = sort(X, 1, 'ascend');
Xrange = max(X)-min(X);
Xprime = X';
U = zeros(S,1);
for c = 1:T-1
   for cc = c+1:T
        d=((Xprime(:, cc) - Xprime(:, c)) - ((cc-c) .* Xrange') ./ T);
        U = U + abs(d);
   end
end
U = T ./ U;

end