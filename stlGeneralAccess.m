function [x, y, z, bdNodes, inNodes] = stlGeneralAccess(filename)
% This function reads an STL file in binary format or ascii format.
% It returns the coordinates of both boundary nodes and interior nodes
% of the 3D polyhedron.
%
% partial REFs:
% 1. stlread() by Doron Harlev for binary STL file via mathswork.com
% 2. https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
%
% Examples:
%        close all
%        clear all
%        clc

%        [x, y, z, bdNodes, inNodes] = stlGeneralAccess('./Region1.stl');
%        if 0
%           patch(x, y, z);
%        else
%           scatter3(bdNodes(:,1), bdNodes(:,2), bdNodes(:,3), 'MarkerFaceColor',[.75 .75 .0]);
%           hold on
%           scatter3(inNodes(:,1), inNodes(:,2), inNodes(:,3),'MarkerFaceColor',[.5 .0 .0]);
%        end
%
% Where
%       filename  ---   Full path of the input STL file.
%       dimSize   ---   [3*1]. The dimSize of coordinate, covering the 3D polyhedron.
%       x         ---   [3*numFacet]. X Position of three points in a triangular facet.
%       y         ---   [3*numFacet]. Y Position of three points in a triangular facet.
%       z         ---   [3*numFacet]. Z Position of three points in a triangular facet.
%       bdNodes   ---   [num_bdNodes*3]. The coordinates of boundary nodes.
%       inNodes   ---   [num_inNodes*3]. The coordinates of interior nodes.
%
% FUNCTION stlGeneralAccess. Version 6.1 Written by JeffZhang. AUG,2020.
% E-mail: jfsufeng[AT]gmail.com || jfzhang2018[AT]zju.edu.cn


if nargin < 1
    error('Too few input arguments')
end
if nargout > 5
    error('Too many output arguments')
end
fid=fopen(filename);
if fid == -1
    error('File could not be opened, pls check the path.')
end
disp('>>>stlGeneralAccess.[START]')
%tic;toc;
nline = 0;
binary = -1;
ascii = -1;
binaryOrText = -1;
stlLines={};
while ~feof(fid)
    tline = fgetl(fid);
    nline = nline + 1;
    if contains(tline, 'endfacet')
        ascii = 1;
        binary = 0;
    end
    if nline >= 50
        break;
    end
end

if ascii == 1 && binary ==0
    disp('FileMode: text');
    binaryOrText = 0;
else
    disp('FileMode: binary');
    binaryOrText = 1;
end

num_facet = 0;
frewind(fid)
switch binaryOrText
    case 0
        while ~feof(fid)
            stlLines=[stlLines;fgets(fid)];
        end
        
        %norm = [];
        ver1=[];ver2=[];ver3=[];
        x=[];y=[];z=[];
        nmRow=[];
        %x=zeros(3,num_facet); y=zeros(3,num_facet); z=zeros(3,num_facet);
        for i = 1:size(stlLines,1)
            %tline = fgetl(fid);
            tl = stlLines{i};
            id = strfind(tl,'facet normal');
            if ~isempty(id)
            %if contains(tl,'facet normal')
                nmRow=[nmRow;i];
                num_facet = num_facet+1;
                %tk = tl(id+13:end); %15 or id+13(better&more stable)
                %nm = strsplit(tk,' ');
                %norm = [norm; str2double(nm)]; %% skip
                
                if i+2 <= size(stlLines,1) && i+3 <= size(stlLines,1) && i+4 <= size(stlLines,1)
                    ver1Str = stlLines{i+2};
                    ver2Str = stlLines{i+3};
                    ver3Str = stlLines{i+4};
                    id1 = strfind(ver1Str,'vertex');
                    id2 = strfind(ver2Str,'vertex');
                    id3 = strfind(ver3Str,'vertex');
                    if isempty(id1) || isempty(id2) || isempty(id3)
                    %if ~contains(ver1Str,'vertex')||~contains(ver2Str,'vertex')||~contains(ver3Str,'vertex')
                        error('Exception0x0815 of vertex.')
                    end                              
                    tk1 = ver1Str(id1+7:end);  %11 or id1+7(better&more stable)
                    tk2 = ver2Str(id2+7:end);
                    tk3 = ver3Str(id3+7:end);
                    nm1 = strsplit(tk1,' ');
                    nm2 = strsplit(tk2,' ');
                    nm3 = strsplit(tk3,' ');
                    ver1 = [ver1; str2double(nm1)];
                    ver2 = [ver2; str2double(nm2)];
                    ver3 = [ver3; str2double(nm3)];
                    x=[x;ver1(end,1) ver2(end,1) ver3(end,1)];  %% For patch validation
                    y=[y;ver1(end,2) ver2(end,2) ver3(end,2)];
                    z=[z;ver1(end,3) ver2(end,3) ver3(end,3)];
                end
            end
        end
        x=x';y=y';z=z';
        
    case 1
        
        ftitle=fread(fid,80,'uchar=>schar'); % Read file title
        num_facet=fread(fid,1,'int32'); % Read number of Facets
        
        fprintf('\nTitle: %s\n', char(ftitle'));
        fprintf('Num Facets: %d\n', num_facet);
        x=zeros(3,num_facet); y=zeros(3,num_facet); z=zeros(3,num_facet);
        ver1=zeros(num_facet,3); ver2=zeros(num_facet,3); ver3=zeros(num_facet,3);
        for i=1:num_facet
            norm=fread(fid,3,'float32'); % normal coordinates, ignored for now
            ver1(i,:)=fread(fid,3,'float32'); % vertex 1
            ver2(i,:)=fread(fid,3,'float32'); % vertex 2
            ver3(i,:)=fread(fid,3,'float32'); % vertex 3
            col=fread(fid,1,'uint16'); % color bytes
            
            x(:,i)=[ver1(i,1); ver2(i,1); ver3(i,1)]; % convert to matlab "patch" compatible format
            y(:,i)=[ver1(i,2); ver2(i,2); ver3(i,2)];
            z(:,i)=[ver1(i,3); ver2(i,3); ver3(i,3)];
        end
        
end
fclose(fid);

bdVertexes=[ver1;ver2;ver3];
bdNodes=unique(bdVertexes, 'rows');

%% imgSize£¬Spacing Infos, demo below. The following interior-node module is additional and optional
%imSize = [512 512 300];
%imSpacing = [1 1 1]; % [.5 0.5 0.5]

minX_bd = min(bdNodes(:,1));
maxX_bd = max(bdNodes(:,1));
minY_bd = min(bdNodes(:,2));
maxY_bd = max(bdNodes(:,2));
minZ_bd = min(bdNodes(:,3));
maxZ_bd = max(bdNodes(:,3));

disp('Default dimSize below.')
dimSize = [ceil(maxX_bd) - floor(minX_bd), ceil(maxY_bd) - floor(minY_bd), ceil(maxZ_bd) - floor(minZ_bd)]  %imSize
imSpacing = (1/10)*dimSize;
for i = 1 : length(imSpacing)
   if imSpacing(i)>1
       imSpacing(i) = 1;
   end
end
imSpacing
origin = [0. 0. 0.];
% posCrd = origin + idxCrd*imSpacing, for medical image computing;
% if using posCrd, stepSize h(i)=imSpacing(i); rangeSz(i) = imSize(i)*imSpacing(i);
inNodesPre=[];
inNodes=[];
for i= floor(minX_bd):1:ceil(maxX_bd)
    for j= floor(minY_bd):1:ceil(maxY_bd)
        for k= floor(minZ_bd):1:ceil(maxZ_bd)
            tmpX = origin(1) + i * imSpacing(1);
            tmpY = origin(2) + j * imSpacing(2);
            tmpZ = origin(3) + k * imSpacing(3);
            if tmpX >= minX_bd && tmpX <= maxX_bd && ...
                    tmpY >= minY_bd && tmpY <= maxY_bd && ...
                    tmpZ >= minZ_bd && tmpZ <= maxZ_bd
                inNodesPre = [inNodesPre; tmpX tmpY tmpZ];
            end
        end
    end
end

if isempty(inNodesPre)
    fprintf('\n!ERROR!Warning:The dimSize or imSpacing may be not suitable.\n')
end

%% fix x, evaluate y,z. (pnpoly method)
for i= 1:1:size(inNodesPre,1)
    polyEdge2dId = find((bdNodes(:,1)>= inNodesPre(i, 1)-.5)&(bdNodes(:,1)<= inNodesPre(i, 1)+.4));
    vertY = bdNodes(polyEdge2dId,2);
    vertZ = bdNodes(polyEdge2dId,3);
    testY = inNodesPre(i,2);
    testZ = inNodesPre(i,3);
    c = 0;
    
    for k=1:1:length(vertY)
        if k == 1
            j = length(vertY);
        else
            j = k - 1;
        end
        if ( ((vertZ(k)>testZ) ~= (vertZ(j)>testZ)) && ...
                (testY < (vertY(j)-vertY(k)) * (testZ-vertZ(k)) / (vertZ(j)-vertZ(k)) + vertY(k)) )
            c = ~c;
        end
    end
    if ~mod(c,2)
        %disp('even -- outside')
    else
        %disp('uneven -- inside')
        inNodes = [inNodes; inNodesPre(i,:)];
    end
end

if isempty(inNodes)
    fprintf('\n!ERROR!Warning:The dimSize or imSpacing, or just the STL file may not be suitable for computing interior nodes.\n')
end

disp('>>>stlGeneralAccess.[OVER]')

end
