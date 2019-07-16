function d=dist_clr(labs1,labs2)
% the Clustering Distance between two label-vectors for the SAME dataset
n=length(labs1);
d=0;
for i=1:n
    for j=1:n
        d=d+ abs( (labs1(i)==labs1(j)) - (labs2(i)==labs2(j)) );
    end
end
d=d/(n^2);
end