%Data analysis of TREQ experimental data for polyA15-CA assembly 


% data = getData;

data = load('CAFinal_2021_02_22.mat');
data = data.data;
T = data{1};
Mc = data{2};
WT = data{3};
WMc = data{4};
C = [7.5,10,12.5,15]*1e-3;
cols = {'r','b','g','m'};

ConstantEnthalpyParameters = [95.4261212282742,0.197526321941171,10.3820826415747,0.322135984791558];
% errors = [0.2,0.0005,0.08,0.02];


IndSitesParameters = [0.507073718403673,-0.00819238833840565,-9.44939618062550,-0.0207135298132286];
% errors = [0.5, .002, 0.05, 0.2];

for i = 1:4
    [T{i},ind] = sort(T{i});
    Mc{i} = Mc{i}(ind);
    WT{i} = WT{i}(ind);
    WMc{i} = WMc{i}(ind);
    
    ind = Mc{i} < log(0.90*25e-6) & Mc{i} > log(0.10*25e-6);
    T{i} = T{i}(ind);
    Mc{i} = Mc{i}(ind);
    WT{i} = WT{i}(ind);
    WMc{i} = WMc{i}(ind);
    
    ind = WT{i} > -1e-5;
    T{i} = T{i}(ind);
    Mc{i} = Mc{i}(ind);
    WT{i} = WT{i}(ind);
    WMc{i} = WMc{i}(ind);
    
    figure(1)
    hold on
    marke = {'<','o','>','square'};
    errorbar(T{i},Mc{i},WMc{i},-WMc{i},WT{i},-1*WT{i},'o','Color',cols{i},'LineWidth',2,'MarkerFaceColor','w','Marker',marke{i})
    plot(T{i},Mc{i},'o','Color',cols{i},'LineWidth',2,'MarkerFaceColor','w','Marker',marke{i},'MarkerSize',8)
end




marke = {'<','o','^','square'};
F1 = figure(1);
           F1.Position = [200 200 600 500];             %positions
           hold on
           hold on
set(gcf,'color','w')
grid off
box on
           ylabel(['ln([M]_c)'])
           xlabel('1/T')
           set(gcf,'color','w')
           set(gca,'LineWidth',2,'FontSize',20,'XColor','k','YColor','k','FontName','Arial')
for i = 1:length(T)
    [t,f] = simDataIndSites(IndSitesParameters,Mc{i},C(i),T{i},WMc{i},WT{i});
    figure(1)
    hold on
   a = plot(t,f,'--','Color',cols{i},'LineWidth',3);
   uistack(a,'bottom');
  [t,f] = simDataConstantEnthalpy(ConstantEnthalpyParameters,Mc{i},C(i),T{i},WMc{i},WT{i});
    figure(1)
    hold on
   a = plot(t,f,'Color',cols{i},'LineWidth',3);
      uistack(a,'bottom');
end

for i = 1:length(T)
    a = plot(T{i},Mc{i},'o','Color',cols{i},'LineWidth',2,'MarkerFaceColor','w','Marker',marke{i},'MarkerSize',8);
    uistack(a,'bottom');
end
for i = 1:length(T)
    a = plot(T{i},Mc{i},'o','Color',cols{i},'LineWidth',2,'MarkerFaceColor','w','Marker',marke{i},'MarkerSize',8);
end

eT = {};
eMc = {};
eW = {};
IndSitesBoot = [];
allT = [];
allMc = [];
allWT = [];
allWMc = [];
inds = 0;

parfor j = 1:500
    eT = {};
eMc = {};
eW = {};
eWT = {};
eWMc = {};
allT = [];
allMc = [];
allWT = [];
allWMc = [];
inds = 0;
for i = 1:4
    allT = [allT;T{i}];
    allMc = [allMc;Mc{i}];
    allWT = [allWT;WT{i}];
    allWMc = [allWMc;WMc{i}];
    inds = [inds;length(allT)];
end
        [~,ind] = datasample(allT,length(allT));
        ind = sort(ind);
        allT = allT(ind);
        allMc = allMc(ind);
        allWT = allWT(ind);
        allWMc = allWMc(ind);
        cc = 1;
        for i = 2:length(inds)
            cind = ind<inds(i) & ind>sum(inds(i-1)+1);
        eT{i-1} = allT(cind);
        eMc{i-1} = allMc(cind);
        eWT{i-1} = allWT(cind);
        eWMc{i-1} = allWMc(cind);
        end
%     [p,diff] = fminsearch(@meritConstantEnthalpy,ConstantEnthalpyParameters,[],eT,eMc,C,eWMc,eWT);
    [p,diff] = fminsearch(@meritIndSites,IndSitesParameters,[],eT,eMc,C,eWMc,eWT);
    IndSitesBoot(j,:) = p;
    j
    
end



eT = {};
eMc = {};
eW = {};
ConstantEnthalpyBoot = [];
allT = [];
allMc = [];
allWT = [];
allWMc = [];
inds = 0;

parfor j = 1:500
    eT = {};
eMc = {};
eW = {};
eWT = {};
eWMc = {};
allT = [];
allMc = [];
allWT = [];
allWMc = [];
inds = 0;
for i = 1:4
    allT = [allT;T{i}];
    allMc = [allMc;Mc{i}];
    allWT = [allWT;WT{i}];
    allWMc = [allWMc;WMc{i}];
    inds = [inds;length(allT)];
end
        [~,ind] = datasample(allT,length(allT));
        ind = sort(ind);
        allT = allT(ind);
        allMc = allMc(ind);
        allWT = allWT(ind);
        allWMc = allWMc(ind);
        cc = 1;
        for i = 2:length(inds)
            cind = ind<inds(i) & ind>sum(inds(i-1)+1);
        eT{i-1} = allT(cind);
        eMc{i-1} = allMc(cind);
        eWT{i-1} = allWT(cind);
        eWMc{i-1} = allWMc(cind);
        end
    [p,diff] = fminsearch(@meritConstantEnthalpy,ConstantEnthalpyParameters,[],eT,eMc,C,eWMc,eWT);
    ConstantEnthalpyBoot(j,:) = p;
    j
    
end




'Constant Enthalpy'
ConstantEnthalpyParams = mean(ConstantEnthalpyBoot)
ConstantEnthalpyErrrors = std(ConstantEnthalpyBoot)

'Independant Sites'
IndSitesParams = mean(IndSitesBoot)
IndSitesErrrors = std(IndSitesBoot)


function rss = meritConstantEnthalpy(params,T,K,C,WMc,WT)
rss = 0;
for i = 1:length(T)
    [t,f] = simDataConstantEnthalpy(params,K{i},C(i),T{i},WMc{i},WT{i});
    rssX = ((T{i}-t(:))./WT{i}).^2;
    rssY = ((K{i}-f(:))./WMc{i}).^2;
    rss = rss+sum((rssX+rssY).^(1/2));
end
% rss
end
function rss = meritIndSites(params,T,K,C,WMc,WT)
rss = 0;
for i = 1:length(T)
    [t,f] = simDataIndSites(params,K{i},C(i),T{i},WMc{i},WT{i});
    rssX = ((T{i}-t(:))./WT{i}).^2;
    rssY = ((K{i}-f(:))./WMc{i}).^2;
    rss = rss+sum((rssX+rssY).^(1/2));
end
% rss
end
function [t,f] = simDataIndSites(params,currK,currC,currT,WMc,WT)
pI = 1;   
H0 = params(pI);
pI = pI+1;
S0 = params(pI);
pI = pI+1;
H1 = params(pI);
pI = pI+1;
S1 = params(pI);
pI = pI+1;


Tref = 25+273.15;
R = 8.314e-3/4.184;
TT = (0:0.01:60)+273.15;

   
Y = exp(-H0./(R.*TT)+S0/R);
K1 = exp(-H1./(R.*TT)+S1/R);
KK = 1./(Y.*(1+K1.*currC).^15);
   f = [];
   t = [];
    for l = 1:length(currK)
%         WT(l) = 1e-11;
    [vc,cind] = min((((currK(l)-log(KK))/WMc(l)).^2+((currT(l)-1./TT)/WT(l)).^2).^(1/2));

    t(l) = 1./TT(cind);
    f(l) = log(KK(cind));
%     figure(1)
% hold on
% plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'g','LineWidth',3)
% 
%     [vc,cind] = min((((1./TT-currT(l))./WT(l)).^2));
% 
%     plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'k')
%     
%         [vc,cind] = min((((log(KK)-currK(l))./WMc(l)).^2));
% plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'r')

    end
%     plot(1./TT,log(KK),'b')
end
function [t,f] = simDataConstantEnthalpy(params,currK,currC,currT,WMc,WT)
pI = 1;   
H0 = params(pI);
pI = pI+1;
S0 = params(pI);
pI = pI+1;
n = params(pI);
pI = pI+1;
Cp = params(pI);
% Cp = 0;
Tref = 25+273.15;
R = 8.314e-3/4.184;
TT = (0:0.01:60)+273.15;
H = H0+Cp*(TT-Tref);
S = S0+Cp*log(TT./Tref);
   KK = currC.^-n.*exp(-(H-TT.*S)./(R*TT));
   f = [];
   t = [];
    for l = 1:length(currK)
%         WT(l) = 1e-11;
    [vc,cind] = min((((currK(l)-log(KK))/WMc(l)).^2+((currT(l)-1./TT)/WT(l)).^2).^(1/2));

    t(l) = 1./TT(cind);
    f(l) = log(KK(cind));
%     figure(1)
% hold on
% plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'g','LineWidth',3)
% 
%     [vc,cind] = min((((1./TT-currT(l))./WT(l)).^2));
% 
%     plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'k')
%     
%         [vc,cind] = min((((log(KK)-currK(l))./WMc(l)).^2));
% plot(currT(l),currK(l),'or')
% % plot(1./TT(cind),log(KK(cind)),'ok')
% % plot(1./TT,log(KK),'b')
% plot([currT(l),1./TT(cind)],[currK(l),log(KK(cind))],'r')

    end
%     plot(1./TT,log(KK),'b')
end
function data = getData()
E = getExperiments();
T = {[],[],[],[]};
Mc = {[],[],[],[]};
WT = {[],[],[],[]};
WMc = {[],[],[],[]};
for j =2:length(E)-1
    c = E{j};
xx = [];
yy = [];
UlineX = c{2}.ux;
UlineY = c{2}.uy;
FlineX = c{2}.lx;
FlineY = c{2}.ly;
for i = 2:length(c)-1
yy = [];
xx = [];
    for jj = 1:1e3
        [UlY,e1] = datasample(UlineY,length(UlineY));
        [FlY,e2] = datasample(FlineY,length(FlineY));
        UlX = UlineX(e1);
        FlX = FlineX(e2);
        p1 = polyfit(FlX,FlY,1);
        p2 = polyfit(UlX,UlY,1);
c{i}.MF = p1(1); 
c{i}.BF = p1(2);
c{i}.MUF = p2(1);
c{i}.BUF = p2(2);
c{i} = getFUF(c{i});
y = c{i}.FUF(:)*25e-6;
x = c{i}.x(:)+273.15;
if c{i}.rate < 0
[m,ind] = max(y);
ind = y>mean(y);
else
  [m,ind] = min(y);  
  ind = y<mean(y);
end
y = y(ind);
x = x(ind);
    [y1,e1] = datasample(y,round(length(y)));
    x1 = x(e1);
    [x1,ind] = sort(x1);
    y1 = y1(ind);
x2 = min(x1):0.01:max(x1);
[p,S,mu] = polyfit(x1,y1,5);
ym = polyval(p,x2,S,mu);
if c{i}.rate < 0
[m,ind] = max(ym);
else
  [m,ind] = min(ym);
end

xx = [xx;(x2(ind))];
yy = [yy;(ym(ind))];
    end
X = 1./mean(xx);
eX = 1/mean(xx) - (1/(mean(xx)-std(xx)));
Y = log(mean(yy));
eY = (log(mean(yy))-log(mean(yy)-std(yy)));
ind = ceil(j/5);
T{ind} = [T{ind};X];
Mc{ind} = [Mc{ind};Y];
WT{ind} = [WT{ind};eX];
WMc{ind} = [WMc{ind};eY];

end
data = {T,Mc,WT,WMc};
end
end
function c = getExperiments()
i = 1;
c = {};
c{i} = autoloadExperiments(classExperiment,'./Data/705mM-1');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/705mM-2');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/705mM-3');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/705mM-4');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/705mM-5');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/10mM-1');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/10mM-2');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/10mM-3');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/10mM-4');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/10mM-5');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/1205mM-1');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/1205mM-2');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/1205mM-3');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/1205mM-4');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/1205mM-5');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/15mM-1');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/15mM-2');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/15mM-3');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/15mM-4');
i = i+1;
c{i} = autoloadExperiments(classExperiment,'./Data/15mM-5');
i = i+1;
end