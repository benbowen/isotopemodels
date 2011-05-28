
function [P,monoiso]=isotope(formula,params,var_el,var_mass,var_abund)
% 
% Implementation is specifically designed in 2010 & 2011 
% by Ben Bowen and Curt Fischer to enable chemical formulae to utilize
% arbitrary nuclide definitions and the entire periodic table.
ele = params.nucleiinfo.ele;
mass = params.nucleiinfo.mass;
abund = params.nucleiinfo.abund;

c = struct2cell(formula);
c(cellfun(@isempty,c)) = {0};
M = cell2mat(c');
% element's isotopes and their natural abundance (Ref:
% http://physics.nist.gov/PhysRefData/Compositions/index.html)
%load isotopic variables 'mass', 'abund' and 'el', containing
%isotope masses, abundances, and elemental symbols for each element



%replace abundances with new abundances
if ~isempty(var_el)
    for i = 1:length(var_el)
        row_n = strmatch(var_el{i},ele,'exact');
        if ~isempty(row_n)
            abund(row_n,:)=0;
            abund(row_n,1:length(var_abund{i})) = var_abund{i};
            if ~isempty(var_mass{i})
                mass(row_n,:)=0;
                mass(row_n,1:length(var_mass{i}))=var_mass{i};
            end
        else
            ele=cat(2,ele,var_el{i});
            temp=zeros(1,size(mass,2));
            mass=cat(1,mass,temp);
            abund=cat(1,abund,temp);
            mass(end,1:length(var_mass{i}))=var_mass{i};
            abund(end,1:length(var_abund{i}))=var_abund{i};
        end
    end
end
syminM = seqmatch(fieldnames(formula),ele,'exact',true);


% element's isotopes and their natural abundance (Ref:
% http://physics.nist.gov/PhysRefData/Compositions/index.html)

%% 
A = mass(syminM,:);
B = abund(syminM,:)/100;
[m mx]=max(B,[],2);
xx=sub2ind(size(A),[1:size(A,1)]',mx);
monoiso=M*A(xx);
% A(:,mx)
% monoiso=0;
% for i = 1:size(A,1)
%     monoiso=monoiso+A(i,mx(i))*M(i);
% end
%% taken from isotop.m which is freely distributed at http://www.ms-utils.org/isotop.html
% 
% calculate isotopic distributions of molecules using the FFT
%
% (c) Magnus Palmblad, 1999
%

MAX_MASS=2^13;      % fast radix-2 fast-Fourier transform algorithm is used
MAX_ELEMENTS=size(A,1);
A2=zeros(MAX_ELEMENTS,MAX_MASS);                 % isotopic abundancies stored in A


for i = 1:size(A,1)
    t=B(i,B(i,:)~=0);
    A2(i,round(A(i,1))+1:round(A(i,1))+length(t))=t;
end


tA=fft(A2,[],2);                     % FFT along each element's isotopic distribution

% this is the slowest part
ptA=ones(1,MAX_MASS);
for i=1:MAX_ELEMENTS,
  ptA=ptA.*(tA(i,:).^M(i));         % multiply transforms (elementwise)
end
% end slowest part

riptA=real(ifft(ptA));              % inverse FFT to get convolutions

id=zeros(1,MAX_MASS);
id(1:MAX_MASS-1)=riptA(2:MAX_MASS); % shift to real mass

% bar(id);   % bar plot of isotopic distributionind=find(ptA>CUTOFF)';
xx=find(id>params.cutoff);
P=[xx' id(xx)'];
% P=id(id>params.cutoff);
