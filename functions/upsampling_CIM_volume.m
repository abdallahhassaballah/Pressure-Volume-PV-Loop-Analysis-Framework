function VV_= upsampling_CIM_volume(VV,n)
    VV_ = nan(n,size(VV,2));
    for ii=1:size(VV,2)
            v = VV(:,ii)';
            x = 1:1:length(v);
            xq = linspace(1,length(v),n);
            Volume_new = interp1(x,v,xq);
            VV_(:,ii) = Volume_new';
    end
end