function terms = basisTerms(x,remDeg,res,terms)

if remDeg==0 || isempty(x)
  terms = [terms,res];
else
  for i=0:remDeg
    new_res = res * x(1)^i;
    terms = basisTerms(x(2:end),remDeg-i,new_res,terms);
  end
end
end