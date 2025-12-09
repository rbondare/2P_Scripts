function [dFF,Ca_baseline]=estimate_continuous_dff(TimeCa,Ca_F,Ca_F_neuropil,neuropil_factor,window,mediantype,varargin)



window_base=ceil(window(1)/mean(diff(TimeCa(1,:))));
window_base=-window_base:window_base;
windowsize_base=numel(window_base);

window_F=ceil(window(end)/mean(diff(TimeCa(1,:))));
window_F=-window_F:window_F;
windowsize_F=numel(window_F);

Ca_Fc=Ca_F-Ca_F_neuropil*neuropil_factor;
Ca_baseline=zeros(size(Ca_F),'single');

switch varargin{1}
    case 'min'
        Ca_baseline=movmin(Ca_Fc,windowsize_base,2,omitnan);
    case 'mean'
        Ca_baseline=movmean(Ca_Fc,windowsize_base,2,omitnan);
    case 'prctile'
        for t=1:size(Ca_F,2)
            cwin=window_base+t;
            cwin(cwin<1 | cwin>size(Ca_F,2))=[];
            Ca_baseline(:,t)=prctile(Ca_Fc(:,cwin),varargin{2},2);
        end
    otherwise
        error('unknown windowing function')
end

dFFtemp=Ca_Fc-Ca_baseline;

switch mediantype
    
    case 'global'
        dFF=dFFtemp./median(dFFtemp,2)-1;
    case 'windowed'
        dFF=dFFtemp./movmedian(dFFtemp,windowsize_F,2)-1;
    otherwise
        error('unknown median type')
end

end
