function  k = gray2int(g)
%  k = gray2int(g)  solves  g = int2gray(k)  for  k ,
%  but only if nonnegative integer  g <= bitmax .  The
%  time taken is proportional to the number of bits
%  needed to hold  g .  See also  grays.m,  int2gray.m,
%  and  graystep.m .             W. Kahan,  8 July 2007

%if(any(any( (g<0)|(g~=round(g))|(g>bitmax) ))),  G = g,
if(any(any( (g<0)|(g~=round(g))|(g>flintmax) ))),  G = g,
    error(' gray2int(G)  needs an array G of small nonnegative integers.')
  end
k = g ;
while any(any(g))
    g = fix(0.5*g) ;  k = bitxor(k, g) ;  end
