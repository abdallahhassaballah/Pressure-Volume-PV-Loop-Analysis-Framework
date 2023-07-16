function cim_sumary = Extract_cim_Summary(data,cim)
    ed_frame = cim.ed_frame;
    nframes  = cim.nframes;
    vol = cim.Volumes;
    
    if ed_frame==1 ||ed_frame==nframes 
        Volume = vol(1:nframes);
    else
        Volume = vol(ed_frame:nframes);
    end
    
    
    cim_sumary = data;
    for KK = 1:size(data,2)    
        PV_lenght = size(data(KK).LVP,1);
        cim_sumary(KK).model_name=cim.model_name;
        cim_sumary(KK).ed_frame = cim.ed_frame;
        cim_sumary(KK).es_frame = cim.es_frame;
        cim_sumary(KK).nframes = cim.nframes;
        cim_sumary(KK).EDV = cim.EDV;
        cim_sumary(KK).ESV = cim.ESV;
        cim_sumary(KK).SV = cim.SV;
        cim_sumary(KK).MASS = cim.MASS;
        cim_sumary(KK).EF = cim.EF;
        cim_sumary(KK).GLS = cim.GLS;
        cim_sumary(KK).HR = cim.HR;
        cim_sumary(KK).Volumes_org = cim.Volumes;
        cim_sumary(KK).Volumes = upsampling_CIM_volume(cim.Volumes', PV_lenght);
    end
end