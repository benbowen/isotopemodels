clear all
close all
clc
%% set required parameters
params.resolution=0.1;
params.windowsize=50;
params.cutoff=1e-4;
params = set_nuclei_composition(params);
% make a molecule with 4 deuterium positions where the percent incorporation is 9% at each position
formula1=[]
% C62H102N14O17S1 
% C36H20N2O4
% C16H35NO2
formula1.C=16;
formula1.H=36;
% formula1.Hn=50;
formula1.O=2;
formula1.N=1;
% formula.S=1;
% formula1.P=4;
% formula2=formula1;

% formula2.Hn=50;
% formula2.H=formula1.H-50;
syminM = seqmatch(fieldnames(formula1),params.nucleiinfo.ele,'exact',true);
params.nucleiinfo.ele=params.nucleiinfo.ele(syminM);
params.nucleiinfo.mass=params.nucleiinfo.mass(syminM,:);
params.nucleiinfo.abund=params.nucleiinfo.abund(syminM,:);
clc
tic
% pat1=fsisotope(formula2,params,{'Hn'},{[1.0078 2.00141]},{[100-0.25 0.25]});
% pat1=fsisotope(formula2,params,{'C'},{[]},{[5 95]});

pat1=fsisotope(formula1,params,{'C'},{[]},{[100-15.2 15.2]});
toc
pat2=fsisotope(formula1,params,[],{[]},{[]});

toc
pat1(:,2)=pat1(:,2)/max(pat1(:,2));
pat2(:,2)=pat2(:,2)/max(pat2(:,2));
stem(pat1(:,1),pat1(:,2),'r.')
hold on
stem(pat2(:,1)+0.01,pat2(:,2),'k.')
hold off
% 916 is where it should be
%% compare two isotopic distributions
pat1(:,1)=pat1(:,1)+randn(size(pat1(:,1)))*0.01;
mztol=0.25;
minnum=2;
F = CompareTwoIsotopomerDistributions(pat1,pat2,mztol,minnum)
%% optimize the isotopic composition of a single element to match a pattern
params.minnum=3;
element='C';
LB=0;
UB=99;
X0=2;
charge=1;
params.plot=1;
params.massrange=7;
x=patternsearch(@(x)OptimizeIsotopicComposition(x,pat1,formula1,charge,params,element),X0,[],[],[],[],LB,UB)
%% make a molecule based on the natural isotopic abundance of elements
params.resolution=0.0001;
params.windowsize=50;
params.cutoff=1e-3;
params = set_nuclei_composition(params);

formula2.C=25;
formula2.H=84;
formula2.O=15;
formula2.N=12;
formula2.P=4;

tic
syminM = seqmatch(fieldnames(formula2),params.nucleiinfo.ele,'exact',true);
params.nucleiinfo.ele=params.nucleiinfo.ele(syminM);
params.nucleiinfo.mass=params.nucleiinfo.mass(syminM,:);
params.nucleiinfo.abund=params.nucleiinfo.abund(syminM,:);

tic

pat1=isotope(formula1,params,{'C'},{[]},{[5 95]});
toc
pat2=isotope(formula1,params,[],{[]},{[]});

toc
pat1(:,2)=pat1(:,2)/max(pat1(:,2));
pat2(:,2)=pat2(:,2)/max(pat2(:,2));
stem(pat1(:,1),pat1(:,2),'r.')
hold on
stem(pat2(:,1),pat2(:,2),'k.')
hold off
%%

[pat2,m]=isotope(formula2,params,[],{[]},{[]});
toc
pat2(:,2)=pat2(:,2)/max(pat2(:,2));
stem(pat2(:,1),pat2(:,2),'r.');
m
%%
pat2=GetMultiNomIsotopicPatternFromFormula(formula2,params,[]);
pat2(:,2)=pat2(:,2)/max(pat2(:,2));
%% plot them
stem(pat1(:,1),pat1(:,2),'r.');
hold on
stem(pat2(:,1),pat2(:,2),'k.');
hold off
%% write an optimization function that calculates n and f
edges=[0:1:10];
X=zeros(size(edges));
for i = 1:length(edges)
    [X(i) F(i)] = fminbnd(@(x)estimateDpattern(x,pat1(:,2),formula2,edges(i),params),0,99.99);
end
%% write an optimization function that calculates n and f
edges=[0:1:10];
X=zeros(size(edges));
for i = 1:10
    temp=pat1(:,2);
    r=randn(size(temp));
    r=r.*temp*0.05;
    temp=temp+r;
    [X(i) F(i)] = fminbnd(@(x)estimateDpattern(x,temp,formula2,4,params),0,99.99)
end
% the differences "pat1(:,2)-pat2(:,2)" shows that the sum of counts in the
% first channel (ie: the monoisotopic channel) has the same counts as the
% rest of the channels.  The implies that the "lost" monoisotopic channel
% is captured in the remainder of the spectra.
%
% the original spectra and convolute it with "something"  The question is 
% a) what pattern and 
% b) what relative weight
%%
close all
formula3=[]
formula3.Hn=1;
pat3=GetMultiNomIsotopicPatternFromFormula(formula3,params,{'Hn'},{[0.01 99.99]});
pat3(:,2)=pat3(:,2)/max(pat3(:,2));
pat3(1,2)=0;
formula4=[]
formula4.H=1;
pat4=GetMultiNomIsotopicPatternFromFormula(formula4,params,[]);
pat4(:,2)=pat4(:,2)/max(pat4(:,2));

% r=pat3(:,2)-pat4(:,2);
% pat3(1,2)=max(pat3(:,2));
% stem(pat4(:,1),r)
%%
% pat3(:,2)=pat3(:,2)/max(pat3(:,2));
% pat4(:,2)=pat4(:,2)/max(pat4(:,2));


%%
close all
% temp=pat2(:,2)-pat4(:,2)+pat3(:,2);
% temp=conv(pat2(:,2),r,'same')
temp=conv(pat2(:,2),pat3(:,2),'same')
% temp=conv(pat2(:,2),pat3(:,2)-pat4(:,2),'same')
% temp=conv(pat2(:,2),pat3(:,2)-pat4(:,2),'same');
temp=temp/max(temp);
stem(pat2(:,1),temp,'g.','markersize',40)
hold on
stem(pat2(:,1),pat1(:,2),'.k');
hold on
stem(pat2(:,1),pat2(:,2),'.r')
hold off