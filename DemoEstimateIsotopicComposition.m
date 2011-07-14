clear all
close all
clc
%%
addpath(genpath('/Users/bpb/matlabtoolbox/'))

%% set required parameters
params.resolution=0.1;
params.windowsize=5;
params.cutoff=1e-4;
params = set_nuclei_composition(params);
formula1=[]
% C36H21N2O4 Scytonemin
formula1.C=36;
formula1.H=21;
formula1.O=4;
formula1.N=2;
% trip unneccessary elements from library
syminM = seqmatch(fieldnames(formula1),params.nucleiinfo.ele,'exact',true);
params.nucleiinfo.ele=params.nucleiinfo.ele(syminM);
params.nucleiinfo.mass=params.nucleiinfo.mass(syminM,:);
params.nucleiinfo.abund=params.nucleiinfo.abund(syminM,:);
% calculate isotopologue distributions
tic
pat1=isotope(formula1,params,{'C'},{[]},{[100-2 2]}); %labeled
toc
tic
pat2=isotope(formula1,params,[],{[]},{[]}); %unlabeled
toc
pat1(:,2)=pat1(:,2)/max(pat1(:,2));
pat2(:,2)=pat2(:,2)/max(pat2(:,2));
stem(pat1(:,1),pat1(:,2),'r.')
hold on
stem(pat2(:,1),pat2(:,2),'k.')
hold off
%% optimize the isotopic composition of a single element to match a pattern
element='C';
LB=0;
UB=99;
X0=1;
charge=1;
params.minnum=5;%problematic for heavily labeled molecules
params.plot=0; %set to 0 to speed up 10x (and not show the plot)
% params.massrange=7;%disabled: problematic for heavily labeled molecules
pat1(:,2)=pat1(:,2)+0.0025*(randn(size(pat1(:,2))).*pat1(:,2));%add white noise
tic
x=patternsearch(@(x)OptimizeIsotopicComposition(x,pat1,formula1,charge,params,element),X0,[],[],[],[],LB,UB)
toc
%% load some real data
load('/Users/bpb/Data/LCMS/desert-crust_C18pos/CrustIsotopeC18Old.mat')
%%
% 274.275	12m04s	C16H35NO2
% 229.2166@1378 is our target C12H27N3O1
formula1=[];
formula1.C=12;
formula1.H=27;
formula1.N=3;
formula1.O=1;
tic
pat1=isotope(formula1,params,{'C'},{[]},{[100-5.2 5.2]}); %labeled
toc
%%
for iii = 1:length(Data)
mz=[Data(iii).features(:).centroidmz];
rt=[Data(iii).features(:).centroidrt];
intensity=[Data(iii).features(:).area];
targetmz=229.2166;
targetrt=1378;
xx1=find(mz>(targetmz-0.1) & mz<(targetmz+3.1) & ...
    (rem(mz,1)-rem(targetmz,1))<0.02 & ...
    abs(rt-targetrt)<30);
measured=zeros(1e6,3);
c=1;
for i = 1:length(xx1)
    L=length(Data(iii).features(xx1(i)).data.mz);
measured(c:c+L-1,1)=round(Data(iii).features(xx1(i)).data.mz);
measured(c:c+L-1,2)=Data(iii).features(xx1(i)).data.rt;
measured(c:c+L-1,3)=Data(iii).features(xx1(i)).data.intensity;
c=c+L;
end
measured(c:end,:)=[];
% measured(measured(:,3)<100,:)=[];
% optimize the isotopic composition of a single element to match a pattern
element='C';
LB=0.1;
UB=10;
X0=1;
charge=1;
params.plot=0; %set to 0 to speed up 10x (and not show the plot)
% params.massrange=7;%disabled: problematic for heavily labeled molecules
tic
u=unique(measured(:,2));
savex{iii}=zeros(length(u),3);
for i = 1:length(u)
    xx=find(measured(:,2)==u(i));
    if length(unique(measured(xx,1)))>1
params.minnum=length(unique(measured(xx,1)));%problematic for heavily labeled molecules
% [x,f]=patternsearch(@(x)OptimizeIsotopicComposition(x,measured(xx,[1 3]),formula1,charge,params,element),X0,[],[],[],[],LB,UB);
[x,f]=fminbnd(@(x)OptimizeIsotopicComposition(x,measured(xx,[1 3]),formula1,charge,params,element),LB,UB);

savex{iii}(i,1)=u(i);
savex{iii}(i,2)=x;
savex{iii}(i,3)=f;
[iii i]
    end
end
end
%%
close all
for i = 1:9
    ax(i)=subplot(3,3,i), plot(savex{i}(:,1),savex{i}(:,2),'.');
    ylim([0 2])
    grid on
end
linkaxes(ax,'xy')
%%
close all
c=[1 0 0;1 0 0;1 0 0;0 1 0;0 1 0;0 1 0;0 0 1;0 0 1;0 0 1]
for i = 1:9
    [s sx]=sort(savex{i}(:,1));
    p=plot(s,smooth(savex{i}(sx,2)),'.');
    set(p,'color',c(i,:))
    hold on
    ylim([0 2])
    grid on
end
% linkaxes(ax,'xy')
%% find some metabolites
mz=[Data(4).features.centroidmz];
rt=[Data(4).features.centroidrt];
intensity=[Data(4).features.height];
idx=[mz<400 & rt>300 & intensity>10000];
scatter(mz(idx),rt(idx),[],log10(intensity(idx)),'filled')

%%
