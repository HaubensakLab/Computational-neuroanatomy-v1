% Copyright (C) 2016 VRVis.
% All rights reserved.
% Contact: VRVis Forschungs-GmbH (office@vrvis.at)
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the VRVis Forschungs-GmbH.
% 4. Neither the name of the VRVis Forschungs-GmbH nor the
%    names of its contributors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%

% This script requires the human connectome project workbench
% http://www.humanconnectome.org/software/connectome-workbench.html
% and adding its libraries to LD_LIBRARY_PATH
% export LD_LIBRARY_PATH=/home/users/ganglberger/workbench-linux64-v1.1.1/workbench/libs_linux64

% This script maps the task fMRI activation maps from HCP_S900_787_tfMRI_ALLTASKS_level3_zstat1_hp200_s2_MSMAll.nii
% (900 Subjects Data Release Reference https://www.humanconnectome.org/study/hcp-young-adult/document/900-subjects-data-release)
% to biopsy-site level from the Allen Human Brain Atlas. Points on the cortex are mapped to the closest biopsy-site
% within 10 mm. This is required since the biopsy sites are not on the surface of the cortex!
% download cifti open file from https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ

hcp_workbench_command = '/home/users/ganglberger/workbench-linux64-v1.1.1/workbench/exe_linux64/wb_command'

tfmriFilename = '..//storage//HCP_S900_787_tfMRI_ALLTASKS_level3_zstat1_hp200_s2_MSMAll.dscalar.nii';

%download cifti open file from https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ
cii = ciftiopen(tfmriFilename,hcp_workbench_command); 

[status,cmdout] = system([hcp_workbench_command ' -file-information ' tfmriFilename]);
cmdout



[status,cmdout] = system([hcp_workbench_command ' -cifti-export-dense-mapping ' tfmriFilename ' COLUMN -volume-all ..//storage//voxelIndizes_TFMRI.txt']);
[status,cmdout] = system([hcp_workbench_command ' -cifti-export-dense-mapping ' tfmriFilename ' COLUMN -surface CORTEX_LEFT ..//storage//verticesIndizesLeft_TFMRI.txt']);
[status,cmdout] = system([hcp_workbench_command ' -cifti-export-dense-mapping ' tfmriFilename ' COLUMN -surface CORTEX_RIGHT ..//storage//verticesIndizesRight_TFMRI.txt']);


cortexLeftFile = 'S900.L.white_MSMAll.32k_fs_LR.surf.gii';
cortexRightFile = 'S900.R.white_MSMAll.32k_fs_LR.surf.gii';


[status,cmdout] = system([hcp_workbench_command ' -surface-coordinates-to-metric ' cortexLeftFile ' ..//storage//coordinatesVertices_left_TFMRI.gii']);
[status,cmdout] = system([hcp_workbench_command ' -surface-coordinates-to-metric ' cortexRightFile ' ..//storage//coordinatesVertices_right_TFMRI.gii']);
posToMNItransform = [-2,0,0,90;0,2,0,-126;0,0,2,-72;0,0,0,1];


 fid = fopen('..//storage//voxelIndizes_TFMRI.txt');
 indizesFile = textscan(fid,'%f%f%f%f','delimiter',' ');
 fclose(fid);
 
 voxelCoordinatesMNI = (posToMNItransform * [[indizesFile{2},indizesFile{3},indizesFile{4}]'; ones(1, size(indizesFile{1},1))])';
 voxelCoordinatesMNI=voxelCoordinatesMNI(1:end,1:3);
 voxelIndizes=indizesFile{1}+1;

 verticesIndizes = [];
 verticesCoordinatesMNI = [];
 
 cFile = gifti('..//storage//coordinatesVertices_left_TFMRI.gii');
 
 fid = fopen('..//storage//verticesIndizesLeft_TFMRI.txt');
 indizesFile = textscan(fid,'%f%f','delimiter',' ');
 fclose(fid);
 verticesIndizes = [verticesIndizes;indizesFile{1}+1];
 verticesCoordinatesMNI = [verticesCoordinatesMNI;cFile.cdata(indizesFile{2}+1,:)];

 cFile = gifti('..//storage//coordinatesVertices_right_TFMRI.gii');
 
 fid = fopen('..//storage//verticesIndizesRight_TFMRI.txt');
 indizesFile = textscan(fid,'%f%f','delimiter',' ');
 fclose(fid);
 verticesIndizes = [verticesIndizes;indizesFile{1}+1];
 verticesCoordinatesMNI = [verticesCoordinatesMNI;cFile.cdata(indizesFile{2}+1,:)];
 

load('..//storage//atlasRegionsHuman.mat');

coordinates = zeros(size(indexX,1),3);
coordinates(:,1)=indexX;
coordinates(:,2)=indexY;
coordinates(:,3)=indexZ;

 
niftiIndizes = [];
niftiIndizesCoordinateIDs = [];

notFoundIndex=[];


coordinatesSurroundings = {};

countIndizesFromVoxel = 0;
countIndizesFromVertices = 0;
couldntFindIndex = 0;
for cI = 1:size(coordinates,1)
    [d,p] = min((voxelCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (voxelCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (voxelCoordinatesMNI(:,3)-coordinates(cI,3)).^2);

    if(sqrt(d)<=(5))
       countIndizesFromVoxel=countIndizesFromVoxel+1;
       niftiIndizes = [niftiIndizes,voxelIndizes(p)];
       niftiIndizesCoordinateIDs = [niftiIndizesCoordinateIDs,(cI)];
       
       coordinatesSurroundings{cI} = [voxelIndizes(sqrt((voxelCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (voxelCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (voxelCoordinatesMNI(:,3)-coordinates(cI,3)).^2)<=5);verticesIndizes(sqrt((verticesCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (verticesCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (verticesCoordinatesMNI(:,3)-coordinates(cI,3)).^2)<=5)];
       
      else
        [d,p] = min((verticesCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (verticesCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (verticesCoordinatesMNI(:,3)-coordinates(cI,3)).^2);
        [d,p] = min((verticesCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (verticesCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (verticesCoordinatesMNI(:,3)-coordinates(cI,3)).^2);
        if(sqrt(d)<=(10))   %%10mm distance necessary since biopsy sites are deeper then the surface vertices!
            countIndizesFromVertices = countIndizesFromVertices+1;
            niftiIndizes = [niftiIndizes,verticesIndizes(p)];
            niftiIndizesCoordinateIDs = [niftiIndizesCoordinateIDs,(cI)];
            
            coordinatesSurroundings{cI} = [verticesIndizes(sqrt((verticesCoordinatesMNI(:,1)-coordinates(cI,1)).^2 + (verticesCoordinatesMNI(:,2)-coordinates(cI,2)).^2 + (verticesCoordinatesMNI(:,3)-coordinates(cI,3)).^2)<=10)];
      
          else
            couldntFindIndex = couldntFindIndex+1;
            notFoundIndex = [notFoundIndex,cI];
            end
    end
end

 disp(sprintf('...got region information vox=#%d / vert=#%d / notF=#%d',countIndizesFromVoxel,countIndizesFromVertices,couldntFindIndex))
 disp(sprintf('...got coordinates information pos=#%d / con row_cols=#%d / notF_Voxels=#%d / notF_Vertices=#%d',length(niftiIndizesCoordinateIDs),length(unique(niftiIndizes)), length(setdiff(voxelIndizes,niftiIndizes)), length(setdiff(verticesIndizes,niftiIndizes))))

 usedContrasts=[81:83,31:33,63:65,37:49,75:77,69:71,1:11,15:22];
 taskfMRI = zeros(size(coordinates,1),length(usedContrasts));

for actContrastIndex = 1:length(usedContrasts)
    taskfMRI(niftiIndizesCoordinateIDs,actContrastIndex)=cii.cdata(niftiIndizes,usedContrasts(actContrastIndex));
end

 

[status,cmdout] = system([hcp_workbench_command ' -file-information ' tfmriFilename ' -only-map-names >> ..//storage//tasknames.txt']);

 fid = fopen('..//storage//tasknames.txt');
 tasknames = textscan(fid,'%s','delimiter','\n');
 fclose(fid);

 tasknames=tasknames{1}(usedContrasts);
 
save('..//storage//taskfMRI.mat','taskfMRI');
save('..//storage//taskfMRI.mat','tasknames','-append');



