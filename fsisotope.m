
function P=fsisotope(formula,params,var_el,var_mass,var_abund)
% Calculation of isotopomer abundance according to
% Rockwood, A. L., Van Orden, S. L. and Smith, R. D. (1996), 
% Ultrahigh Resolution Isotope Distribution Calculations. 
% Rapid Communications in Mass Spectrometry, 10: 54?59. 
% doi: 10.1002/(SICI)1097-0231(19960115)10:1<54::AID-RCM444>3.0.CO;2-Z
% 
% Implementation is specifically designed in 2010 & 2011 
% by Ben Bowen and Curt Fischer to enable chemical formulae to utilize
% arbitrary nuclide definitions and the entire periodic table.

c = struct2cell(formula);
c(cellfun(@isempty,c)) = {0};
M = cell2mat(c');
M=cat(2,M,0);
syminM = seqmatch(fieldnames(formula),params.nucleiinfo.ele,'exact',true);
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


c = struct2cell(formula);
c(cellfun(@isempty,c)) = {0};
M = cell2mat(c');
M=cat(2,M,0);
syminM = seqmatch(fieldnames(formula),ele,'exact',true);
% element's isotopes and their natural abundance (Ref:
% http://physics.nist.gov/PhysRefData/Compositions/index.html)

%% 
A = mass(syminM,:);
B = abund(syminM,:)/100;
A=cat(1,A,zeros(1,size(A,2)));

A(end,1)=1;
atomicnumber=round(A(:,1))-1;
xx=find(A==0);
A=A-atomicnumber*ones(1,size(A,2));
A(xx)=0;
atomicnumber(end)=-1;

B=cat(1,B,zeros(1,size(B,2)));
B(end,1)=1;
A=round(A*1e5);
A=cat(3,A,B);
%% taken from isotop_fs.m which is freely distributed on 
% Calculates isotopic distributions including isotopic fine structure
% of molecules using FFT and various scaling 'tricks'.
%
% (C) 1999 by Magnus Palmblad, Division of Ion Physics, Uppsala Univ.
% Acknowledgements:
% Lars Larsson-Cohn, Dept. of Mathematical Statistics, Uppsala Univ.,
% for help on theory of convolutions and FFT.
% Jan Axelsson, Div. of Ion Physics, Uppsala Univ. for comments and ideas
% Contact Magnus Palmblad at magnus.palmblad@gmail.com if you should
% have any questions or comments.

MAX_ELEMENTS=length(M);
MAX_ISOTOPES=size(A,2);
if isfield(params,'cutoff')
CUTOFF = params.cutoff;
else
CUTOFF=1e-10;       % relative intensity cutoff for plotting
end

WINDOW_SIZE=params.windowsize; %input('Window size (in Da) ---> ');

RESOLUTION=params.resolution;%input('Resolution (in Da) ----> ');  % mass unit used in vectors
if RESOLUTION < 0.00001  % minimal mass step allowed
    RESOLUTION = 0.00001;
elseif RESOLUTION > 0.5  % maximal mass step allowed
    RESOLUTION = 0.5;
end

R=0.00001/RESOLUTION;  % R is used to scale nuclide masses (see below)

temp=squeeze(A(:,1,1))*R;
temp(end)=0;
Mmi=round(temp)'*M';

WINDOW_SIZE=WINDOW_SIZE/RESOLUTION;   % convert window size to new mass units
WINDOW_SIZE=2^nextpow2(WINDOW_SIZE);  % fast radix-2 fast-Fourier transform algorithm

if WINDOW_SIZE < round(496708*R)+1
    WINDOW_SIZE = 2^nextpow2(round(496708*R)+1);  % just to make sure window is big enough
end

% mass shift so Mmi is in left limit of window:

FOLDED=floor(Mmi/(WINDOW_SIZE-1))+1;  % folded FOLDED times (always one folding due to shift below)

% shift distribution to 1 Da from lower window limit:
M(MAX_ELEMENTS)=ceil(((WINDOW_SIZE-1)-mod(Mmi,WINDOW_SIZE-1)+round(100000*R))*RESOLUTION);
MASS_REMOVED=atomicnumber'*M';  % correction for 'virtual' elements and mass shift

ptA=ones(WINDOW_SIZE,1);
for i=1:MAX_ELEMENTS
    tA=zeros(WINDOW_SIZE,1);
    for j=1:MAX_ISOTOPES
        if A(i,j,1) ~= 0
            tA(round(A(i,j,1)*R)+1)=A(i,j,2);  % put isotopic distribution in tA
        end
    end
    
    tA=fft(tA) ;  % FFT along elements isotopic distribution  O(nlogn)
    
    tA=tA.^M(i);  % O(n)
    
    ptA=ptA.*tA;  % O(n)
end

ptA=real(ifft(ptA));  % O(nlogn)

MA=linspace((FOLDED*(WINDOW_SIZE-1)+1)*RESOLUTION+MASS_REMOVED,(FOLDED+1)*(WINDOW_SIZE-1)*RESOLUTION+MASS_REMOVED, WINDOW_SIZE-1)';

ind=find(ptA>CUTOFF)';
P=[MA(ind) ptA(ind)];

