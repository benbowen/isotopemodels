clear all
close all
clc
%% set required parameters
params.resolution=0.001;
params.windowsize=4;
params = set_nuclei_composition(params);
% make a molecule with 4 deuterium positions where the percent incorporation is 9% at each position
formula1=[]
formula1.C=25;
formula1.H=80;
formula1.Hn=4;
formula1.O=15;
formula1.N=12;
formula1.P=4;
clc
tic
pat1=fsisotope(formula1,params,{'Hn','C'},{[],[0.5 1 1.5]},{[100-9 9],[90 7.5 2.5]});
toc
pat1(:,2)=pat1(:,2)/max(pat1(:,2));
stem(pat1(:,1),pat1(:,2))
% 916 is where it should be
%% make a molecule based on the natural isotopic abundance of elements
formula2.C=25;
formula2.H=84;
formula2.O=15;
formula2.N=12;
formula2.P=4;
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