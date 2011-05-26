function biny=GetMultiNomIsotopicPatternFromFormula(formula,params,var_el,var_abund)

% %%
% params.Purity12C=0.011;
% params.Purity13C=0.968;%pe(mx);
% params.Purity15N=0.95;
% params.PercentInc13C=0.975;
% params.PercentInc15N=0.975;
% params.Iso{1}.ID='Carbon';
% params.Iso{1}.p=[0.9890 0.0110];
% params.Iso{1}.m=[0   1.003355];
% params.Iso{2}.ID='Hydrogen';
% params.Iso{2}.p=[0.99985 0.00015];
% params.Iso{2}.m=[0   1.006277];
% params.Iso{3}.ID='Nitrogen';
% params.Iso{3}.p=[0.9963 0.0037];
% params.Iso{3}.m=[0   0.997035];
% params.Iso{4}.ID='Oxygen';
% params.Iso{4}.p=[0.99757 0.00038 0.00205];
% params.Iso{4}.m=[0   1.004216   2.004244];
% params.Iso{5}.ID='Sulfur';
% params.Iso{5}.p=[0.9502 0.0075 0.0421 0.0002];
% params.Iso{5}.m=[0   0.999387   1.995796   3.995007];
% params.Iso{6}.ID='Phosphorus';
% params.Iso{6}.p=[1];
% params.Iso{6}.m=[0];
% %%
% params.H=1.007825; %mass of hydrogen
% params.O=15.994915;
% params.N=14.003074;
% params.S=31.972072;
% params.C=12;
% params.P=30.973763;
% params.ppm_filter=5;
% 
%     massedges=pubchemmass(r(i)):1.003355/5:pubchemmass(r(i))+65;
%     biniso = GetIsotopicPatternFromFormula(pubchemfull(r(i),:),params,massedges,pubchemmass(r(i)));

%This function will crash if mzvec ~= edges(1) within narrow margin
% if abs(edges(1)-min(mzvec))>0.5 | abs(edges(1)-max(mzvec))>0.5
% error('you did not set mzedges to start at mzvec')
% end
c = struct2cell(formula);
c(cellfun(@isempty,c)) = {0};
M = cell2mat(c');
% ws = warning('off','Bioinfo:seqmatch:StringNotFound');
syminM = seqmatch(fieldnames(formula),params.nucleiinfo.ele,'exact',true);
% warning(ws);
% element's isotopes and their natural abundance (Ref:
% http://physics.nist.gov/PhysRefData/Compositions/index.html)
%load isotopic variables 'mass', 'abund' and 'el', containing
%isotope masses, abundances, and elemental symbols for each element

ele = params.nucleiinfo.ele;
mass = params.nucleiinfo.mass;
abund = params.nucleiinfo.abund;

%replace abundances with new abundances
if ~isempty(var_el)
    for i = 1:length(var_el)
    row_n = strmatch(var_el{i},ele,'exact');
    abund(row_n,1:length(var_abund{i})) = var_abund{i};
    end
end

% elements = [mass, abund];
A = mass(syminM,:);
B = abund(syminM,:);
monoisomass=M*A(:,1);
edges=0:params.resolution:params.windowsize;
%%

biny=zeros(length(edges),2);
% for i = 1:size(fmat,1)
    for j = 1:length(M)
        temp=MultinomialIsotopeDist(M(j),B(j,1:4)/100,A(j,1:4)-A(j,1),edges);
        if j==1
            pattern=temp;
        else
            pattern=conv(pattern,temp);
            pattern=pattern(1:numel(temp));
            %             pattern=pattern+temp;
        end
    end
    biny(:,2)=pattern;
    biny(:,1)=edges+monoisomass;
% end

%%
function [y]=MultinomialIsotopeDist(n,p,m,edges)
if n~=0
    x1 = 0:n;
    if length(p)==1
        be=x1(:);
    elseif length(p)==2
        [X1,X2] = ndgrid(x1,x1);
        be=[X1(:) X2(:)];
    elseif length(p)==3
        [X1,X2,X3] = ndgrid(x1,x1,x1);
        be=[X1(:) X2(:) X3(:)];
    elseif length(p)==4
        [X1,X2,X3,X4] = ndgrid(x1,x1,x1,x1);
        be=[X1(:) X2(:) X3(:) X4(:)];
    end
    xx=find(sum(be,2)~=n);
    be(xx,:)=[];
    if length(p)==1
        Y=1;
        x=0;
    else
        Y = mnpdf(be,p);
        x=be*m';
    end
    y=mHist(x,edges-edges(1),Y);
else
    Y=1;
    x=0;
    y=mHist(x,edges-edges(1),Y);
end