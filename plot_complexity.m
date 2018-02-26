Nt=1:10;
Nm=1:11;
P=zeros(length(Nt),length(Nm));
for i=1:length(Nt)
    for j=1:length(Nm)
        [i,j]
        s=0;
        for t=0:min(Nt(i),Nm(j))
            s=s+nchoosek(Nt(i),t)*nchoosek(Nm(j),t)*factorial(t);
        end
        
        P(i,j)=s;
        
    end
end


figure
bar3(P)

set(gca,'ZScale','log')

llim = .1;

h = get(gca,'Children');

for i = 1:length(h)
    ZData = get(h(i), 'ZData');
    ZData(ZData==0) = llim;
    set(h(i), 'ZData', ZData);
end

% set(gca,'XTickLabel',{'GLgn3','CUT4-U'})
% view(119,24)
xlabel('#Targets')
ylabel('#Measurements')
% zlabel('Number of points')
plot_prop_paper
