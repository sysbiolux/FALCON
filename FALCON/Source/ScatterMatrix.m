function [AllCorrelations] = ScatterMatrix(Dataset, Names, Category, Legend, Plot)


mn=size(Dataset,2);
Dens=log(size(Dataset, 1)/max(Category));

D=Dataset;
AllCorrelations=nan(mn,mn,max(Category));

colors=distinguishable_colors(max(Category));
leg=1;
if Plot
    figure,
end

for c=1:mn
    for cc=c:mn
        if Plot
            subplot(mn, mn, (c-1)*mn + cc)

            if c==cc
                hist(D(:,c),40)
%                 xlim([0 1]);
            end
        end
        if ~(c==cc)
            for cat=1:max(Category)
                DataX=D(Category==cat,c);
                DataY=D(Category==cat,cc);
                Fit=fitlm(DataX, DataY);
                Rsq=Fit.Rsquared.Adjusted;
                AllCorrelations(c,cc,cat)=Rsq;
                if Plot
                    scatter1=scatter(DataX, DataY,[],colors(cat,:), '.'); hold on,
                    if numel(DataX)>100
                        scatter1.MarkerFaceAlpha = 1/Dens;
                        scatter1.MarkerEdgeAlpha = 1/Dens;
                    end
                end

            end
        end
        if Plot
%             ylim([0 1]);
%             xlim([0 1]);
        end
%         Fit=fitlm(D(:,c), D(:,cc));
%         Rsq=Fit.Rsquared.Adjusted;
        if Plot
            if ~c==cc

                text(0.5,0.5,num2str((ceil(Rsq*100))/100))
            end
        end

        if Plot
            if c==1 && cc==mn && leg==1
                leg=0;
                legend(Legend,'Location','eastoutside')
            end
        end

        if Plot
            if c==1
                title(Names{cc})
                set(gca, 'FontSize', 5);
            end

            if cc==1
                title(Names{c})
                set(gca, 'FontSize', 5);
            end
        end
        if Plot
            set(gca,'xticklabels',[])
            set(gca,'yticklabels',[])
        end
        if ~exist('Rsq'), Rsq=NaN; end
        if isnan(Rsq), Rsq=0; end
        disp([num2str(c),'/', num2str(mn),'   ',  num2str(cc),'  ', repmat('x',1,ceil(Rsq*10))])
    end
    if Plot
        drawnow
    end
end
% if Plot
%     set(gcf,'Position',get(0,'Screensize'));
% end

end