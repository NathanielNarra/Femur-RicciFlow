function [b1,b2] = Detect_Boundaries(TR)
% A slightly modified version of code from 'iso2mesh' package
% Credit to original author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>

V = TR.Points;
F = TR.ConnectivityList;
ed = freeBoundary(TR);

loops=[];
loops=[loops,ed(1,:)];
loophead=ed(1,1);
loopend=ed(1,end);
ed(1,:)=[];

while(length(ed))
    idx=[find(ed(:,1)==loopend)',find(ed(:,2)==loopend)'];
%    if(length(idx)>1) error('self intersecting curve is unsupported'); end
    if(length(idx)==1)
        idx=idx(1);
        newend=setdiff(ed(idx,:),loopend);
        if(newend==loophead)
            loops=[loops,nan];
            ed(idx,:)=[];
            if(size(ed,1)==0) break; end
            loops=[loops,ed(1,:)];
            loophead=ed(1,1);
            loopend=ed(1,end);
            ed(1,:)=[];
            continue;
        else
            loops=[loops,newend];
        end
        loopend=newend;
        ed(idx,:)=[];
    end
end

t = bwlabel(~isnan(loops));
b1 = (loops(t==1))';
b2 = (loops(t==2))';
