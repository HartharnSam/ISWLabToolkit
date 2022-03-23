function [im]=dfireadvel(fn,inp)
% DFIREAD - read DigiFlow floating point format (.dfi) images
%
% im=dfiread('filename.dfi'); will read the full DigiFlow image file
% and output the result to the structure IM. This structure contains
% all the "tags" from the DigiFlow file. The image can be found in
% im.cdata and the colormap in im.cmap (256x3 array-RGB).
%
% If you just want to extract the image (assuming only 1 plane exists),
% use
% >> im=dfiread('filename.dfi','simple');
%
% Display the image using IMAGESC since they are 32-bit (which is not
% supported by IMSHOW). To display your image on screen type:
% >> imagesc(im.cdata);
%
% Currently DFIREAD only supports uncompressed images.
%
% Note that images are "flipped" so that their axis' start in the
% upper left corner (in Digiflow they start at the lower left).
%
% For multi-plane images, the im.cdata will be a 3-dimensional matrix
% containing each image plane. The description of each plane can be
% found in the im.Descr cell-array:
% >> im.Descr{1}
% ans =
%
% u
%
% ...which implies that the first plane (im.cdata(:,:,1)) is actually
% the U-component of a vector field. Assuming the second plane is the
% vertical component, we can display the vector field:
% >> quiver(im.cdata(:,:,1),im.cdata(:,:,2));
%
% Alternatively you can display each image plane as an image:
% >> imagesc(im.cdata(:,:,3)), title([strcat(im.Descr{3}),' - ',strcat(im.Contains{3})])
%
% See also: IMAGESC
%
% J.K.Sveen@damtp.cam.ac.uk, Copyright
% Distributed under the GNU GPL license.

% releases:
%
% Version 0.45, Apr 14 -05
% Added a test to avoid errors when encountering unknown tags. This
% may happen when DigiFlow overwrites an existing .dfi file which on
% some DF versions will leave parts of the old file trailing after
% the recently written one.
%
% Version 0.44, Mar 31 -05
% Speedup improvements and numerous small fixes
%
% Version 0.44, Mar 30 -05
% Added small fix to allow plotting if no output is called; in other
% words, if you type
% >> dfiread('myimage.dfi');
% The image will be plotted in a figure
% Also fixed another occurence of uint32 which should be float32
% (caxis statement)
%
% Version 0.43, Nov 19 -04
% Bug fix - some images were read as uint32 when they should have
% been float 32.
%
% Version 0.42, May 25 -04
% The im.Contains field now contains a descripton of what each image
% plane contains (for multi-plane images only). This information
% supplements the info in the im.Descr field. Example
% >> im=dfiread('myimage.dfi');
% >> [strvcat(im.Descr{:}),strvcat(im.Contains{:})]
% This will display the contents. These fields do NOT exist for single
% plane images.
%
% Version 0.41, May 25 -04
% Another update for the multi-plane images. There was an
% undocumented feature in the file-format for the 4108-tag. This had
% no effect on the image reading, though, but meant that we were
% reading 32 bytes of rubbish integers at the end of this tag. This
% is actually a string containing the original colorscheme name.
%
% Version 0.4, May 24 -04
% Swapped caxis elements so that you can do
% >> imagesc(p.cdata); caxis(p.caxis)
% also fixed reading of multi-plane images and fixed the
% image-flipping to also work in this case.
%
%
% Version 0.3, May 10 -04
% corrected 8 and 64 bit reading bug in version 0.1 + small change in
% the id.Format reading + removed output to screen if called from
% overloaded function

% Check input - abort if image does not exist
if ~exist(fn,'file') || exist(fn,'dir')
    error(join(['Could not find "', fn, '"']));
end

if nargin==1, inp='full';
else
    if ~ischar(inp) || ~strcmp(inp,'simple')
        %... - abort if second argument is not 'full' or 'simple'
        disp('Error. Wrong input to DFIREAD. Please check your input.')
        return
    end
end

%check if this is called from an overloaded function or prompt.
%Suppress output if this the former...
[sta,~]=dbstack;

% Open file
f = fopen(fn,'r');

% Read the first tags
p.idFormat=fread(f,32,'*char')'; %char(fread(f,32,'uchar'))';
p.Version=(fread(f,1,'uint32'));

DataTypes={'1001','cdata','uint8' ;
    '11001','cdata','uint8';
    '1004','cdata','float32';
    '11004','cdata','float32';
    '1008','cdata','float64';
    '11008','cdata','float64';
    '1014','caxis','float32';
    '1018','caxis','float64';
    '1100','rescale','uint32';
    '12001','cdata','uint8';
    '12004','cdata','uint8';
    '12008','cdata','uint8';
    '2000','cmap','uint8';
    '2001','cmapname','uint8';
    '2002','cmapnamevar','uchar';
    '3000','Descr','uchar'; ...
    '3001','Comment','uchar';
    '3002','CreatingProcess','uchar';
    '3018','ImageTime','float32';
    '4008','ImageCoordinates','uint32';
    '4108','ImagePlaneDetails','uint32'};

containstag={'000','None';
    '001','Grayscale';
    '002','Red';
    '003','Green';
    '004','Blue';
    '10','Scalar field';
    '101','xCoordinate';
    '102','yCoordinate';
    '103','zCoordinate';
    '201','xVector';
    '202','yVector';
    '203','zVector';
    '301','Unknown?';
    '302','Unknown?';
    '303','Unknown?';
    '311','dxVector/dx';
    '312','dxVector/dy';
    '313','dxVector/dz';
    '321','dyVector/dx';
    '322','dyVector/dy';
    '323','dyVector/dz';
    '331','dzVector/dx';
    '332','dzVector/dy';
    '333','dzVector/dz';
    };


%TAG is made up from datatype and nbytes
p.DataType=(fread(f,1,'uint32'));

%TAG is followed by data of type DataType.
% We now read TAG followed by data until DataType is empty, i.e. no
% more data in the file.
while ~isempty(p.DataType) %isstr(dec2hex(p.DataType)),
    
    p.nBytes=(fread(f,1,'uint32'));
    
    %consider trying to replace strvcat with a cell-array
    [indy,indx]=find(strcmp(char(dec2hex(p.DataType)),DataTypes));
    if ~isempty(indy) && ~isempty(indx)
        str=DataTypes{indy,2};
        typ=DataTypes{indy,3};
        swstr=DataTypes{indy,1};
    else
        str='UNKNOWN';
        swstr='UNKNOWN';
    end
    if size(sta,1)<=1
        %disp(['Reading datatype: #',dec2hex(p.DataType),' - ',str])
    end
    switch lower(swstr)
        case {'12001','12004','12008'}
            % Compressed images are not supported.
            if size(sta,1)<=1
                disp(['Error. DFIREAD does not support compressed images. Please', ...
                    ' reload your image(s) in DigiFlow and edit stream to be', ...
                    ' saved without compression.']); return % hence we abort here
            end
            
            p.nx=fread(f,1,'uint32');
            p.ny=fread(f,1,'uint32');
            p.nz=fread(f,1,'uint32');
            %myimage=zeros(p.nx,p.ny,p.nz);
            p.szCompress=fread(f,1,'uint32');
            eval(['p.c=(fread(f,p.szCompress,''',typ,'''))'';']);
            %How about writing this to a .mat file manually since that uses
            %gzip compression? hmmm....that would be slow
            
            %None of the attempts below work off the shelf
            %save -ascii tmp.dat.zip
            %!unzip tmp.dat.zip
            %p.c=load('tmp.dat','ascii');
            %test2=num2cell(p.c);
            %test3=Huff06(test2);
            %eval(['p.',str,'=mxuncompress(myimage(:)'',length(myimage(:)''),p.c(:)'',length(p.c(:)''));'])
            %eval(['p.',str,'=reshape(p.',str,',p.ny,p.nx,p.nz);'])
            
        case {'1001','1004','1008'}
            % here we read the image array - 2D image type
            p.nx = fread(f,1,'uint32');
            p.ny = fread(f,1,'uint32');
            eval(['p.',str,' = (fread(f,(p.nBytes-8)/4,''',typ,'''))'';']);
            eval(['p.',str,' = reshape(p.',str,',p.ny,p.nx);'])
        case {'11001','11004','11008'}
            % here we read the image array - 3D image type (although the
            % third dimension may be == 1
            if strcmp(dec2hex(p.DataType),'11001'), divi = 1;
            elseif strcmp(dec2hex(p.DataType),'11004'), divi = 4;
            else, divi = 8;
            end
            
            p.nx = fread(f,1,'uint32');
            p.ny = fread(f,1,'uint32');
            p.nz = fread(f,1,'uint32');
            eval(['p.',str,' = (fread(f,(p.nBytes-12)/divi,''',typ,'''))'';']);
            eval(['p.',str,' = (reshape(p.',str,',p.nx,p.ny,p.nz));']);
            
        case {'1014','1018'}
            % coloraxis data etc
            cax = fread(f,2,typ);
            p.caxis(1:2) = cax;%(2); %num2str(cax(1:p.nBytes/2));
            %p.caxis(2)=cax(1); %num2str(cax(1+p.nBytes/2:end));
            
        case {'1100'}
            %image rescaling stuff - ignore for now in matlab
            p.nxWant=fread(f,1,'uint32');
            p.nyWant=fread(f,1,'uint32');
            p.method=fread(f,1,'uint32');
            
        case {'2000'}
            % Colormap
            test = fread(f,768,'uint8');
            test = test(:);
            p.cmap(:,1)=test(1:256); %fread(f,256,'uint8'); % Red component
            p.cmap(:,2)=test(257:512); %fread(f,256,'uint8'); % Green component
            p.cmap(:,3)=test(513:768); %fread(f,256,'uint8'); % Blue component
            
        case {'2001','2002'}
            % Colormap-name, ignore for now in Matlab
            if strmp(str,'2001')
                iLen = 64;
            else
                iLen = fread(f.iLen,'uint32');
            end
            p.cmapname = fread(f,iLen,'*char')';%char(fread(f,iLen,'uchar'))';
            
        case {'3000','3001','3002'}
            %...and here char-arrays, comments etc
            eval(['p.',str,'=char(fread(f,p.nBytes,''',typ,'''))'';']);
            
        case {'3018'}
            % Timing stuff
            p.iFrame=fread(f,1,'uint32');
            p.Reserved=fread(f,1,'uint32');
            p.Time=fread(f,1,'double');
            p.tStep=fread(f,1,'double');
            p.tFirst=fread(f,1,'double');
            
        case {'4008'}
            % world coordinate stuff
            p.Kind = fread(f,1,'uint32');
            p.xWorldPerPixel = fread(f,1,'double');
            p.yWorldPerPixel = fread(f,1,'double');
            p.xOriginWorld = fread(f,1,'double');
            p.yOriginWorld = fread(f,1,'double');
            p.xUnits = fread(f,16,'char')';%char(fread(f,16,'uchar'))';
            p.yUnits = fread(f,16,'char')';%char(fread(f,16,'uchar'))';
            p.OriginalName = fread(f,64,'*char')';%char(fread(f,64,'uchar'))';
            
        case {'4108'}
            % 3D image plane stuff
            p.nPlanes=fread(f,1,'uint32');
            for planes=1:p.nPlanes
                p.Contains{planes}=dec2hex(fread(f,1,'uint32'));
                p.PlaneDescr{planes}=(fread(f,32,'char=>char'))';
                p.ParamA(planes)=fread(f,1,'double');
                p.ParamB(planes)=fread(f,1,'double');
                p.ParamC(planes)=fread(f,1,'double');
                p.ParamD(planes)=fread(f,1,'double');
                p.OrigColorSchemeName=(fread(f,32,'char=>char'))';
                %p.ParamE(planes)=fread(f,1,'double');
                %p.ParamF(planes)=fread(f,1,'double');
                %p.ParamG(planes)=fread(f,1,'double');
                %p.ParamH(planes)=fread(f,1,'double');
                %p.Contains1=fread(f,1,'uint32');
                %p.Descr1=char(fread(f,32,'uchar'))';
            end
            %    case {'UNKNOWN'}
            %     % not much to do in this case but to continue.
            %      dummybytes=fread(f,1,'uint32');
            %      dummybytes=fread(f,dummybytes,'uint32');
        otherwise
            %disp('Error'); %return
            msgstring='There were some unknown tags in your .dfi file.';
    end
    
    p.DataType=(fread(f,1,'uint32'));
    
end

fclose(f);

% Now clean up the structure
p = rmfield(p,'nBytes');
p = rmfield(p,'DataType');

% ...and rotate/flip image to correspond to Matlab/Image-processing
% standard (axis ij)
[ix,iy,iz] = size(p.cdata);
myim = nan(iy, ix, iz);
if iz>1
    for i=1:iz
        myim(:,:,i)=flipud(p.cdata(:,:,i)');
    end
else
    myim=flipud(p.cdata');
end
p.cdata=[];
p.cdata=myim;

% Subsequently change the p.Contains from hexadecimal values to
% understandable text, using the conversion-table CONTAINSTAG
if isfield(p,'Contains')
    for i=1:length(p.Contains)
        % consider changing strvcat to using a cell array
        coind = strcmp(containstag,char(p.Contains{i}));
        p.Contains{i} = containstag{coind,2};
    end
    % generate coordinate system - this will be slow for interpolated
    % vector fields
end
if isfield(p, 'xOriginWorld')
    %imsize = [p.nx*p.ParamB(1)  p.ny*p.ParamB(2)];
    X = linspace(p.xOriginWorld, (p.xOriginWorld + p.xWorldPerPixel*(p.nx-1)), p.nx);
    Y = linspace(p.yOriginWorld, (p.yOriginWorld + p.yWorldPerPixel*(p.ny-1)), p.ny)';
    %X = p.ParamB(1):p.ParamC(1):imsize(1);
    %Y = p.ParamB(3):p.ParamC(3):imsize(2);
    X=repmat(X,p.ny,1);
    Y=repmat(flipud(Y(:)),1,p.nx);
    %convert to world coordinates
    %X=p.xWorldPerPixel*X + p.xOriginWorld;
    %Y=p.yWorldPerPixel*Y + p.yOriginWorld;
    p.x=X;
    p.y=Y;
end
if strcmp(inp,'simple')
    if size(sta,1)<=1 % suppress output if this is called from within
        % another file. Else display to screen
        disp('Truncating structure.... simple output')
    end
    im=p.cdata;
else
    im=p;
end

if nargout==0
    % Here we draw the image in a figure window if no output is
    % specified. We need to treat single plane images, velocity-fields
    % and multi-plane images separatly.
    disp('Opening figure to display contents....')
    if size(p.cdata,3)==1
        imagesc(p.cdata), caxis([p.caxis])
        %title([strcat(im.Descr{1}),' - ',strcat(im.Contains{1})])
    elseif size(p.cdata,3)==2
        % this should be a vector field
        quiver(1:10:p.nx, 1:10:p.ny, p.cdata(1:10:end,1:10:end,1), ...
            p.cdata(1:10:end,1:10:end,2),'k')
        title(['Vectors show (',strcat(im.Descr{1}),',',strcat(im.Descr{2})])
    elseif size(p.cdata,3)==7
        % This should be a PTV or PIV file (?) - we'll only plot vectors with
        % vorticity as background
        imagesc(p.cdata(:,:,3)), caxis([p.caxis])
        title(['Vectors show (',strcat(im.Descr{1}),',',strcat(im.Descr{2}),...
            '). Background is ',strcat(im.Descr{3})])
        hold on
        quiver(1:10:p.nx, 1:10:p.ny, p.cdata(1:10:end,1:10:end,1), ...
            p.cdata(1:10:end,1:10:end,2),'k')
        colorbar
        hold off
    end
end


% finally warn if we encountered any unknown tags
%if exist('msgstring','var')
%     disp(['Warning! ',msgstring]);
%end