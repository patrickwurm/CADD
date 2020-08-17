function delta = kronecker( dimension )
delta = tenzeros( dimension );
for i=1:dimension
  for j=1:dimension
    if( i == j )
      delta(i,j)=1;
    end
  end
end
end