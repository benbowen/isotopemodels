function F = OptimizeIsotopicComposition(x,pat1,formula,charge,params,element)
pat2=fsisotope(formula,params,{element},{[]},{[100-x x]});
F = CompareIsotopomerDistributionToTheoretical(pat2,pat1,charge,params); %pass in theoeretical first
x