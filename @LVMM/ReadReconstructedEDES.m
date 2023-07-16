function Lmm = ReadReconstructedEDES(Folder, casename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function was added to add recontructed frames after
%Bias correction
%author: Kat Gilbert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lmm = LVMM;
%Read and add ED frame
Files = dir(sprintf('%s/*%s*',Folder, casename));
if size(Files,1) ==2
    Lmm = LVMM;
    Lmm.name = str2num(casename);
    %Lmm.nframes =2;
    for mi = 1:2
        P = ReadCIMModelFileRaw(strcat(Folder,'\',Files(mi).name));
        Lmm.focalLengths(mi) = P.FocalLength;
        Lmm.lambdas(:,mi) = P.LambdaParameters;
        Lmm.mus(:,mi) = P.MuParameters;
        Lmm.thetas(:,mi) = P.ThetaParameters;
        frame = split(Files(mi).name,'_');
        idx = size(frame,1);
        if  contains(frame{idx},'ED')
            Lmm.ed = mi;
        else
            Lmm.es  = mi;
        end
    end
    %Voldata = Lmm.ComputeVolumes();
    %Lmm.EDV = Voldata.endo(Lmm.ed);
    %Lmm.ESV = Voldata.endo(Lmm.es);
    %Lmm.MASS = Voldata.mass;
else
    disp('Error - theres not two files for this case')
    return
end
end