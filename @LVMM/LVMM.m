classdef (ConstructOnLoad = true) LVMM < hgsetget
% LVMM = LV Motion Model, it's a faster version

    % user changeable properties
    properties
        name = '';
        ed = 0;
        es = 0;

        data = struct;                  % this is any data derived from the model
    end
    
    % read-only properties
    properties( SetAccess = private )
        EDV = NaN;            % endocardial volume at ED
        ESV = NaN;            % endocardial volume at ES
        SV = NaN;             % stroke volume
        EF = NaN;             % ejection fraction
        MASS = NaN;           % LV mass at ED
    end
    
    % read-only hidden properties
    properties( SetAccess = private, Hidden = true )
        focalLengths = [];              % 1x(nframes) focal lengths
        lambdas = [];                   % 134x(nframes) lambda parameters
        mus = [];                       % 40x(nframes) mu parameters
        thetas = [];                    % 40x(nframes) theta parameters
        phase_data = struct;            % any data for each phase
    end
    
    % on-demand properties
    properties( Dependent = true, Transient = true )
        nframes;        % number of frames
    end

    methods
        % ------------------- CONSTRUCTOR -------------------
        function L = LVMM(varargin)
        % Class constructor
        %
        %   L = LVMM(model_folder);
        %
        % Author: Avan Suinesiaputra - Center for Advanced Imaging, Univ. of Auckland (2012)
        
        end
        % ------------------- END OF CONSTRUCTOR -------------------
        
        % MAIN FUNCTION: model(idx) --> LV4x4SurfaceModel creation on the fly
        function Li = model(L,i)
            if ( L.nframes == 0 ), Li = LV4x4SurfaceModel;
            elseif( i < 1 || i > L.nframes )
                error('Index out of bound.');
            else
                Li = LV4x4SurfaceModel.Create(L.focalLengths(i),L.lambdas(:,i),L.mus(:,i),L.thetas(:,i));
                if( numel(L.phase_data)>=i )
                    Li.data = L.phase_data(i);
                end
            end
        end
        
        % GET FUNCTIONS: to get on-demand properties
        
        function n = get.nframes(L)
            n = numel(L.focalLengths);
        end
        
        function L = set.ed(L,i)
            L.ed = i;
            L.UpdateGlobalLVFunctions;
        end

        function L = set.es(L,i)
            L.es = i;
            L.UpdateGlobalLVFunctions;
        end
        
        % Movie player
        Play(L,varargin);
        MyPlay(L,varargin);
        PlayAHA(L,varargin);
        
        % Show methods
        h = PlotVolumes(L,varargin);
        
        % Volumes
        V = ComputeVolumes(L);
        V = ComputeRegionalVolumes(L);

        % Other functions
        idx = ShiftEDToFirst(L);
        EstimateEDES(L);
        %ReorderFrames(L,frames);
        Rescale(L,new_focalLengths);
        ExportModel(L,model_folder,varargin);
        
    end

    % these methods do not require an object to be defined first
    methods( Static = true )
        % object creation factory
        L = ReadFromCIMFolder( folder, varargin );
        L = Create( focalLengths, lambdaParams, muParams, thetaParams, varargin );
        L = ReadReconstructedEDES( ModelPath,CurrentCase);
    end
    
    % private methods
    methods( Access = private )
        UpdateGlobalLVFunctions(L);
    end
    
end