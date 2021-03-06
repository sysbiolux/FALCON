function sinaplot(X)

% sinaplot(X) produces a plot of the array X, in such a way that
% every column in the array is represented by a different scatterplot, 
% evidencing the distribution of the datapoints, even when the number of 
% datapoints differs between the groups.
% X is a M-by-P matrix with maximum M samples and P groups
% X can contain NaNs
% uses the excellent distinguishable_colors(P) to choose the right hues
%
% Sebastien  De Landtsheer November 2018
% latest update February 2019
% sebdelandtsheer@gmail.com
DoIt = 1;

if isempty(X)
    X = randn(1000, 3) + [1 2 3];
    disp('No data provided: plotting dummy variables')
end
if ~ismatrix(X) %is this a 2D matrix?
    warning('Too many dimensions');
    DoIt = 0;
end
[m, p] = size(X);

if p < 2 %are there at least two groups to plot?
    warning('There is only one group');
    DoIt = 0;
end
if m < 5 %are there at least 5 values?
    warning('Not enough samples to plot');
    DoIt = 0;
end

MaxDatapoints = max(sum((~isnan(X)))); %The maximum number of non-NaN datapoints in any column
MinDatapoints = min(sum((~isnan(X)))); %The minimum number of non-NaN datapoints in any column
if MinDatapoints < 5
    warning('Not enough samples to plot');
    DoIt = 0;
end

if DoIt
    jit = (rand(size(X))-0.5) * 0.95; %Uniform jitter centered on 0 
    xdef = repmat((1:p), m, 1); %Background X-axis position

    %Figuring out the density of points
    X = sort(X);
    Dens = zeros(m, p);
    k = 2; Dens(k, :) = 3 ./ ((X(k+1, :) - X(k-1, :)));
    k = m-1; Dens(k,:) = 3 ./ ((X(k+1, :) - X(k-1, :)));

    for k = 3:m-2
        Dens(k, :) = 5 ./ ((X(k+2, :) - X(k-2, :)));
    end
    if MaxDatapoints > 16
        for Lim = 4:floor((min(40, sqrt(MaxDatapoints)))/2)
            for k = (Lim+1):(m-Lim)
                Dens(k, :) = (Lim * 2 + 1) ./ ((X(k+Lim, :) - X(k-Lim, :)));
            end
        end
    end


    Dens = Dens ./ max(Dens(:)); %normalizing density
    xval = xdef + (jit .* Dens); %new X-axis values
    Colors = distinguishable_colors(p);
    for Group=1:p
        plot(xval(:, Group), X(:, Group), '.k'), hold on,
        plot([Group-0.3, Group+0.3],[nanmean(X(:, Group)), nanmean(X(:, Group))], '-k', 'Linewidth', 1)
    end
    set(gca, 'XGrid', 'on')
    set(gca, 'XTick', 1:p),
end



