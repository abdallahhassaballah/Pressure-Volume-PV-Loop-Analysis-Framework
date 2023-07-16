classdef (ConstructOnLoad = true) LV4x4SurfaceModel < hgsetget
% A class to generate sampling points for 4x4 LV surface FEM

    % user changeable properties
    properties
        nSamples = [9 9]; % [theta_samples mu_samples]
        nSurfaces = 2;    % define how many surfaces
        
        % this is identity matrix
        T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];  % RC points transformation matrix
        
        % this is for any data that you can attach to a model
        data = struct;
    end
    
    % read-only properties
    properties( SetAccess = private )
        % you can change these model parameter by using SetParams function
        lambdaParams    % the transmural component of the LV model
        focalLength     % focal length of the LV model
        muParams        % the azimuthal component of the LV model
        thetaParams     % the circumferential component of the LV model
        
        G2E             % global to element mapping for each prolate spheroidal components

        surfacePoints   % surface points in rectangular Cartesian
                        % columns are [X Y Z]
        
        surfacePSPoints % surface points in prolate spheroidal coordinates
                        % columns are [theta mu lambda]
                        
        basis           % basis functions

        E   % sample values between element nodes
            % E is N x 4 matrix, where N is the total number of
            % samples and each column is:
            %  - col 1 = sample values for dimension 1 (theta)
            %  - col 2 = sample values for dimension 2 (mu)
            %  - col 3 = surface sample values
            %  - col 4 = element numbers
              
        map   % anatomical mapping of the surface elements
        
        nodes % nodal elements: 4x16, where each columns are the 4 corner of the element
              % in the following order: [bottom-right top-right bottom-left top-left]
    end
    
    % properties that their values are calculated on demand
    properties( Dependent = true )
        nCirc     % number of circumferential sample points
        nAzimuth  % number of azimuthal sample points       
    end
    
    methods
        % ------------------ CONSTRUCTOR -----------------------------------
        function L = LV4x4SurfaceModel(varargin)
        % Class constructor
        %
        %   L = LV4x4SurfaceModel('opt1',val1,'opt2',val2,...);
        %
        % Calling L = LV4x4SurfaceModel will create an empty model with
        % default properties.
        %
        % Available options:
        %   - 'nSamples', [nDim1 nDim2].
        %     Define the number of samples between two element nodes on
        %     the LV surface. Default is [9 9].
        %
        %   - 'nSurfaces', num
        %     Define the number of surfaces (max 2). Default is 2.
        %
        %   - 'initFromCIM', 'cim_file'
        %     Read the model from a CIM model file. See ReadFromCIMModel method.
        %
        %   - 'copyFrom', another_LV4x4SurfaceModel_object
        %     Make a copy from another object.
        %
        % Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)
            
            % get options
            validopts = {'nSamples','nSurfaces'};
            addOpts.initFromCIM = '';
            addOpts.copyFrom = [];
            otherOpts = {};
            for i=1:2:length(varargin)
                if( find(strcmp(varargin{i},validopts)) )
                    L.(varargin{i}) = varargin{i+1};
                elseif( isfield(addOpts,varargin{i}) )
                    addOpts.(varargin{i}) = varargin{i+1};
                else
                    otherOpts = [otherOpts varargin(i:i+1)];
                end
            end
            
            % update samples
            L.UpdateSampling;
            
            % check init option
            if( ~isempty(addOpts.initFromCIM) )

                L.ReadFromCIMModel(addOpts.initFromCIM);
            elseif( ~isempty(addOpts.copyFrom) )
                if( ~isa(addOpts.copyFrom,'LV4x4SurfaceModel') )
                    error('Cannot copy from the given object.');
                end
                L.data = addOpts.copyFrom.data;
                L.G2E = addOpts.copyFrom.G2E;
                L.lambdaParams = addOpts.copyFrom.lambdaParams;
                L.muParams = addOpts.copyFrom.muParams;
                L.thetaParams = addOpts.copyFrom.thetaParams;
                L.focalLength = addOpts.copyFrom.focalLength;
                L.nSamples = addOpts.copyFrom.nSamples;
                L.nSurfaces = addOpts.copyFrom.nSurfaces;
                L.UpdateSampling;
                
            end

        end
        % ------------------ END OF CONSTRUCTOR --------------------------------
       
        % set/get methods
        
        function set.nSamples(L,newNumberOfSamples)
            if( any(newNumberOfSamples<1) ), error('Invalid number of samples.'); end
            if( numel(newNumberOfSamples)~=2 ), error('Dimension of the number of samples must be 2.'); end
            
            L.nSamples = newNumberOfSamples(:)';
            L.UpdateSampling;
        end
        
        function set.nSurfaces(L,newNumberOfSurfaces)
            if( numel(newNumberOfSurfaces) ~= 1), error('Invalid number of surfaces.'); end
            if( newNumberOfSurfaces < 1 ), error('Invalid number of surfaces.'); end
            
            L.nSurfaces = newNumberOfSurfaces;
            L.UpdateSampling; 
        end
        
        function SetParams(L,newFocalLength,newLambdas,newMus,newThetas)
            L.focalLength = newFocalLength;
            L.lambdaParams = newLambdas;
            L.muParams = newMus;
            L.thetaParams = newThetas;
            L.UpdateSampling;
        end

        function n = get.nCirc(L)
            if( isempty(L) ), n = 0;
            else n = size(L.map,2); end
        end
        
        function n = get.nAzimuth(L)
            if( isempty(L) ), n = 0;
            else n = size(L.map,1); end
        end
        
        function P = get.surfacePoints(L)
            % always apply transformation matrix
            P = [L.surfacePoints ones(size(L.surfacePoints,1),1)] * L.T;
        end
                
        % visualization methods in 3D space
        h = plot(L,varargin);
        h = PlotSurfacePoints(L,varargin);
        h = PlotBoundingBox(L,varargin);
        h = PlotWireframe(L,varargin);
        h = PlotSurface(L,varargin);
        
        % visualization methods in element space
        h = PlotElementGrid(L,si,varargin);
        h = PlotElementXi(L,si,P,varargin);
        
        % specialized visualization methods
        h = PlotPolarLambda(L,si,varargin);
        h = PlotPolarRC(L,si,varargin);
        
        % advanced visualization with interaction
        view(L);
        h = PlotSliceVolume(L,V,ax,pos,varargin);
        %BrowseSlices(L);
        
        % input/output functions
        L = ReadFromCIMModel(L,cim_file);
        WriteAsCIMModel(L,cim_file);
        
        % compute points on the surface
        [P,idx] = GetApicalPoints(L);
        [P,idx] = GetBasalPoints(L);
        [P,idx] = GetLongitudinalPoints(L,ci);
        [P,idx] = GetCircumferentialPoints(L,li);
        [P,idx] = GetSurfacePoints(L,si);
        
        % volumetric properties
        BB = GetBoundingBox(L,varargin);
        [V,org,res] = GetVolumeLabel(L,varargin);
        P = GetPointsInsideSurface(L,si,varargin);
        [M,Lx] = GetOuterMask(L,varargin);
        
        % compute surface/triangulation
        F = GetSurfaceFaces(L,si);
        [FV,iTopFaces] = GetClosedSurfaceTriangulation(L,si);
        FV = GetMyocardialTriangulation(L);
        
        % points derived from an LV object
        [P,Fidx] = GetIntersectionWithPlane(L,P0,N0,varargin);
        Ps = ProjectPointsToSurface(L,si,P);
        P = GetIntersectionWithDICOMImage(L,dimg,varargin);
        
        % get the surface element from a given points in cartesian coordinates
        Xi = FindSurfaceElementXi(L,si,P);
        TM = Xi2ThetaMu( L, si, Xi, Ei );
        P = GetSurfacePointsFromElements( L, E);
        
        % global parameters
        [N,idx] = GetNodalParameters(L,si);
        TML = GetGlobalPSPosition(L);
        
        % volume & mass
        V = CalculateVolume(L,si);
        function V = GetEndoVolume(L)
            V = L.CalculateVolume(1);
        end
        function V = GetEpiVolume(L)
            V = L.CalculateVolume(L.nSurfaces);
        end
        mass = GetMass(L,varargin);
        
        % test methods
        [IN,ON] = InsideSurface(L,si,P,varargin);
        
        % REGIONAL ANALYSIS
        % ------------------------------------------------------------
        [mu_levels,idx,P] = GetAHASliceLevels(L);
        [R,mus,thetas] = GetAHAPoints(L,varargin);
        [FV,Pb] = GetAHASurfacePatches(L,G);
        
        [V,FV] = CalculateRegionalAHAVolumes(L,G);
        
        PlotAHASegmentsXi(L,si,varargin);
        % ------------------------------------------------------------
        
        
    end
    
    % these methods do not require an object to be defined first
    methods( Static = true )
        E = RectangularSampling(nElements, nSamples);
        B = ComputeBasisFunctions(elements, nBasis, basisType);
        XYZ = Prolate2Cart( focalLength, lambda, mu, theta );
        TML = Cart2Prolate( focalLength, X, Y, Z );
        
        % object creation factory
        L = Create( focalLength, lambdaParams, muParams, thetaParams, varargin );
        L = CreateEndoEpiSurfacesByTML( focalLength, TML, varargin );
    end

end