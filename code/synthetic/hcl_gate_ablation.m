function T = hcl_gate_ablation(H,D)
% HCL_GATE_ABLATION
% Uses the SAME final evidence, po, frozen mask, and recovered tensor.
% Only the conditional attribution rule changes.
%
% This isolates attribution quality from reconstruction.

E = H.final;
po = E.po;
variants = {'FullHCL','Flat','EntryOnly','NoBlock','NoSlice'};

rows = repmat(struct('Variant','','LevelMacroF1',NaN, ...
                     'EntryF1',NaN,'BlockF1',NaN,'SliceF1',NaN), ...
              numel(variants),1);

for i=1:numel(variants)
    [qe,qb,qs] = shares(E.ae,E.ab,E.as,variants{i});
    piE=po.*qe; piB=po.*qb; piS=po.*qs;
    det=po>=0.50;
    pred=zeros(size(det),'uint8');
    A=cat(4,piE,piB,piS);
    [~,lab]=max(A,[],4);
    pred(det)=uint8(lab(det));

    rows(i).Variant=variants{i};

    if ~isempty(D.truth.level)
        [rows(i).EntryF1,rows(i).BlockF1,rows(i).SliceF1, ...
         rows(i).LevelMacroF1] = level_scores(pred,D.truth.level);
    else
        rows(i).LevelMacroF1=NaN;
    end
end

T=struct2table(rows);
end

function [qe,qb,qs]=shares(ae,ab,as,name)
switch upper(name)
    case 'FULLHCL'
        qs=as;
        qb=(1-as).*ab;
        qe=(1-as).*(1-ab);

    case 'FLAT'
        den=ae+ab+as+eps;
        qe=ae./den; qb=ab./den; qs=as./den;

    case 'ENTRYONLY'
        qe=ones(size(ae)); qb=zeros(size(ae)); qs=zeros(size(ae));

    case 'NOBLOCK'
        qs=as; qb=zeros(size(as)); qe=1-as;

    case 'NOSLICE'
        qs=zeros(size(as)); qb=ab; qe=1-ab;

    otherwise
        error('Unknown gate variant: %s',name);
end
end

function [fe,fb,fs,macro]=level_scores(pred,truth)
fe=class_f1(pred,truth,1);
fb=class_f1(pred,truth,2);
fs=class_f1(pred,truth,3);
present=unique(double(truth(:))); present(present==0)=[];
all=[fe fb fs];
macro=mean(all(present),'omitnan');
end

function f=class_f1(pred,truth,c)
tp=nnz(pred==c & truth==c);
fp=nnz(pred==c & truth~=c);
fn=nnz(pred~=c & truth==c);
if nnz(truth==c)==0
    f=NaN;
else
    f=2*tp/max(2*tp+fp+fn,1);
end
end
