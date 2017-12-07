function [distance]=bp_clr_knn(type,bp_repeat,algo,X,cmin,cmax,k_knn,varargin)
% [cmin,cmax] is the range of cluster number ccc
% distance is a "(cmax-cmin+1) * times" matrix, every number in matrix is
%   the clustering distance between two bootstrap clusters
% when type==1, every bootstrap clusters are compared, so times=repeat*(repeat-1)/2
% when type~=1, every time two new bootstrap samples are drawn, so times=repeat

s=RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

if type==1
    distance=zeros(cmax-cmin+1,bp_repeat*(bp_repeat-1)/2);
    for ccc=cmin:cmax
        bp=cell(1,bp_repeat);
        labs=cell(1,bp_repeat);
        for i=1:bp_repeat
            [bp{i},labs{i}]=bp_onestep(algo,X,ccc,varargin);
        end
        number_now=1;
        for i=1:bp_repeat
            for j=1:(i-1)
                cmn_labs1=knnclassify(X,X(bp{i},:),labs{i},k_knn);
                cmn_labs2=knnclassify(X,X(bp{j},:),labs{j},k_knn);
                distance(ccc-cmin+1,number_now)=dist_clr(cmn_labs1,cmn_labs2);
                number_now=number_now+1;
            end
        end
    end
else
    distance=zeros(cmax-cmin+1,bp_repeat);
    for ccc=cmin:cmax
        for i=1:bp_repeat
            [bp1,labs1]=bp_onestep(algo,X,ccc,varargin);
            [bp2,labs2]=bp_onestep(algo,X,ccc,varargin);
            
            cmn_labs1=knnclassify(X,X(bp1,:),labs1,k_knn);
            cmn_labs2=knnclassify(X,X(bp2,:),labs2,k_knn);
            distance(ccc-cmin+1,i)=dist_clr(cmn_labs1,cmn_labs2);
        end
    end
end
    