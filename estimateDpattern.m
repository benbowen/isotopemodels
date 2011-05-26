function F = estimateDpattern(x,y,formula,hn,params)
formula.H=formula.H-hn;
formula.Hn=hn;
pat=GetMultiNomIsotopicPatternFromFormula(formula,params,{'Hn'},{[100-x x]});
pat(:,2)=pat(:,2)/max(pat(:,2));
F=sum((y-pat(:,2)).^2);
F
x