function [cmn_bp,cmn_labs1,cmn_labs2]=common(bp1,bp2,labs1,labs2)
% intermidiate function
% returns the common bootstrap samples of bp1,bp2
% and the two vectors of label of the common sample
n=length(labs1);
cmn_bp=zeros(n,1);
cmn_labs1=zeros(n,1);
cmn_labs2=zeros(n,1);

i=1; j=1; number=0;
while i<=n && j<=n
    if bp1(i)==bp2(j)
        number=number+1;
        cmn_bp(number)=bp1(i);
        cmn_labs1(number)=labs1(i);
        cmn_labs2(number)=labs2(j);
        i=i+1;
        j=j+1;
    else
        if bp1(i)<bp2(j)
            i=i+1;
        else
            j=j+1;
        end
    end
end
cmn_bp=cmn_bp(1:number);
cmn_labs1=cmn_labs1(1:number);
cmn_labs2=cmn_labs2(1:number);
end
