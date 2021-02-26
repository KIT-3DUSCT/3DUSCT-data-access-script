%%%minimal USCT II data & metadata load and usage script
%% Version 1.4 M. Zapf, KIT 2016-2021
%% Home repository https://github.com/KIT-3DUSCT/3DUSCT-data-access-script
%% Please submit patches, bugs and pull requests :), thank you!
%% License of this code is permissive BSD 3 clause.
%%
%% Optional support libraries: https://github.com/KIT-3DUSCT/SAFT-CPU
%%                             https://github.com/KIT-3DUSCT/Faster-1D-XCorr
%%
%% profile call: "profile off, profile on, USCT_data_convert(); profile off, profile report,"
%% tested with: _exp0066 , _exp119_Pr9_P2 , _exp131_ , exp0010 and all the open sozrce data

function [imag]=USCT_data_convert()
close all
clear all
%Pathdata='/home/katz/work/Hiwi/exp_006_2D_USCT/AScans' %'C:\brustphantom'; %your path to data, please change

%Pathdata='Y:\Data\_USCT3Dv2\Mannheim\exp0010_mannheim_gelatine\brustpute' %Pathdata = '/home/katz/fser/Hiwi/data/Experiments/exp_010_2D_analyticBreast/AScans';
%Pathdata='Y:\Data\_USCT3Dv2\_exp131_measurments_for_neslihan\neslihan_phantom1'
%Pathdata='Y:\Data\_USCT3Dv2\_exp119_Pr9_P2'
Pathdata='Y:\Data\_USCT3Dv2\_exp0066_phantom (MF1.6MHz) (fixed)'

%%% load meta data and constants
load([Pathdata filesep 'info.mat']);
try load(['.' filesep 'geometryFileUSCT3Dv2_3.mat'])
catch
    load([Pathdata filesep 'geometry.mat']) %legacy alternative
end
load([Pathdata filesep 'CE.mat']);
load([Pathdata filesep 'TASTempComp.mat'],'TASTemperature');
load([Pathdata filesep 'Movements.mat'])

%%% experiment flags block - select data subsets
numTAS = 24; %all: 157, 1-157 are possibilities for USCT II
numEmit= 2; %all: 4, 1-4 are possibilities for USCT II
numRec = 3; %all: 9, 1-9 are possibilities for USCT II
useMPs = 1; %experiment dependend: 1-10 typically

%%% signal processing
upsampling=20e6; %%in Hz, upsampling of AScans -> set higher than 10e6, multiplies of 10e6 
matchedFiltering = 1; % pulse compression with the coded excitation signal, increases SNR
eoffset = 20e-7; %in sec, Excitation signal delay, empirical gathered for USCT II
trytosuppressbaddata = 1; %try to detect errornous AScan data via std-dev
angleEmRecComb = 120; %in degree, filters AScans out under specific angles, increases SNR: 180° is everything ("off"), 120° is common USCT transmission suppression

%%%SAFT imageing
envelopeImaging = 1; % 0 = phase coherent imaging, 1 = envelope imaging (more robust against system errors but lower resolution)
SAFT_MEX = 0; % 0 = default SAFT Matlab-code, 1 = faster SAFT MEX-code (needs https://github.com/KIT-3DUSCT/SAFT-CPU), 2 = both for debug
visualization = 1; % 0 = no visualization, 1 = debug visualization
image_save = 0; %stores images on HDD
%%
%%%%end of parameter flags block -> image definition extra block later!

%%%var init
soundvelocity = 1480; %initial value, later updated
SF = 0; %sample frequency, later updated

%%geometry preparation
geomRecBase=[]; geomEmitBase=[]; geomRecBaseNormals=[]; geomEmitBaseNormals=[];
for i=1:numTAS geomRecBase=cat(3,geomRecBase,TASElements(i).receiverPositions);
    %geomRecBaseNormals=cat(3,geomRecBaseNormals,TASElements(i).receiverPositionsNormals);
end
for i=1:numTAS geomEmitBase=cat(3,geomEmitBase,TASElements(i).emitterPositions);
    %geomEmitBaseNormals=cat(3,geomEmitBaseNormals,TASElements(i).emitterPositionsNormals);
end

middlePoint=squeeze(mean(mean(geomRecBase(:,:,:),1),3));
try middlePoint(3)=squeeze(geomRecBase(5,3,37)); %middlePoint(1:2)=[0 rec_pos(37,5,2)];
catch
    middlePoint(3)=mean(mean(geomRecBase(:,3,:),3)); 
end

%%image definition block
imagXYZ=[1024 1024 1];  %X Y Z pixel amount = select any amount in any dimension
imag=zeros(imagXYZ);
imag2=zeros(imagXYZ);
imagsum=zeros(imagXYZ);
imagRes=0.24/imagXYZ(1); %in m (cubes only!)
imagstart=[-0.11 -0.11 0.04]; %[-0.11 -0.11 0.03];%[-0.042 -0.042 0.03]; %in m
imagPosX=imagstart(1):imagRes:imagstart(1)+(imagXYZ(1)-1)*imagRes;
imagPosY=imagstart(2):imagRes:imagstart(2)+(imagXYZ(2)-1)*imagRes;
imagPosZ=imagstart(3):imagRes:imagstart(3)+(imagXYZ(3)-1)*imagRes;

%%for vectorized Matlab SAFT
imagPos=imagstart;
x=imagPos(1):imagRes:imagPos(1)+(imagXYZ(1)-1)*imagRes;
y=imagPos(2):imagRes:imagPos(2)+(imagXYZ(2)-1)*imagRes;
z=imagPos(3):imagRes:imagPos(3)+(imagXYZ(3)-1)*imagRes;
imagPosAll=cat(4,repmat(x',[1 imagXYZ(2) imagXYZ(3)]), repmat(y,[imagXYZ(1) 1 imagXYZ(3)]),repmat(shiftdim(z,-1),[imagXYZ(1) imagXYZ(2) 1]));

%%downsampled image mesh
downsampling=10; %ratio
[meshx,meshy,meshz]=meshgrid(imagPosX(1:downsampling:end),imagPosY(1:downsampling:end),imagPosZ(1:downsampling:end));

if SAFT_MEX==1 addsig2vol_3(8); end %%init 8 threads (set for your CPU number) in SAFT MEX

i=0;

%loop over all data
for Mp=1:min(size(MovementsListreal,1),useMPs)
    movement=MovementsListreal(Mp,:);
    rotshift=0; %shift by 20 degree ->
    if any(isnan(movement))
        transform_matrix = eye(4);       
    else
        transform_matrix = makehgtform('zrotate',2*pi*(rotshift+movement(1))/360)*makehgtform('translate',[0 0 movement(2)]);
    end
    
    %transform it to the actual lift and rot-pos
    geomRec=zeros(size(geomRecBase)); geomEmit=zeros(size(geomEmitBase));
    for i=1:size(geomRecBase,1)
        for j=1:size(geomRecBase,3)
            temp=([geomRecBase(i,:,j) 1]) * transform_matrix';
            geomRec(i,:,j)=[temp(1)/temp(4) temp(2)/temp(4) temp(3)/temp(4)];
            %temp=([geomRecBaseNormals(i,:,j) 1]) * transform_matrix';
            %geomRecBaseNormals(i,:,j)=[temp(1)/temp(4) temp(2)/temp(4)];
        end
    end
    for i=1:size(geomEmitBase,1)
        for j=1:size(geomEmitBase,3)
            temp=([geomEmitBase(i,:,j) 1]) * transform_matrix';
            geomEmit(i,:,j)=[temp(1)/temp(4) temp(2)/temp(4) temp(3)/temp(4)];
            %temp=([geomEmitBaseNormals(i,:,j) 1]) * transform_matrix';
            %geomEmitBaseNormals(i,:,j)=[temp(1)/temp(4) temp(2)/temp(4)];
        end
    end
    
    for eT=1:numTAS
        for eE=1:numEmit
            clear AScanDatatype
            %%data reconstruction
            %%[Gain,Data]=loadAscan_v2(eT,eE,rT,rE,Mp,Pathdata); %load single Data slow 
            load(sprintf('%s%sTAS%03d%sTASRotation%02d%sEmitter%02d.mat',Pathdata,filesep,eT,filesep,Mp,filesep,eE));
                     
            
            try AScanDatatype, catch, AScanDatatype='int16'; end %if not exist assume old legacy measurement with int16 (PS: "try catch" is faster than "exist")
            
            try load([Pathdata filesep 'CEMeasured.mat']); %measured Coded excitation
                CE=CEMeasured;
            catch,
                load([Pathdata filesep 'CE.mat']); %defined input coded exciation CE CS_SF 
            end
            
            if strcmp(AScanDatatype,'float16')
                AScans=convertfp16tofloat(AScans);
                if exist('CEMeasured','var') %reconstruct only IF MEASURED
                   CE=double(convertfp16tofloat(CE));
                end
            end
            %flip dim of CE to same dim as AScans
            if size(CE,1)==1 && length(CE)>1 CE=CE'; end
                        
            if Bandpassundersampling==1                
                if exist('CEMeasured','var') %reconstruct only IF MEASURED
                    if size(CE,1)<size(AScans,1) CE(size(AScans,1),:)=0; end %%padding for reconstructBPS
                    CE=ReconstructBandpasssubsampling(double(CE));                    
                end  
                AScans=ReconstructBandpasssubsampling(double(AScans)); %now 10MHz
                %now 10MHz
                SF=10e6;
            else
                %already 10MHz 
                AScans=double(AScans);
                SF=10e6;
            end
            
            %%data check -> try to suppress bad data blocks
            if trytosuppressbaddata == 1
                val=std(double(AScans)./repmat(max(AScans),[size(AScans,1) 1]));
                
                %set Ascans to 0 which have std==0 , constant
                ACV=0;
                AScans(:,val==ACV) = 0;
                
                %set Ascans to 0 which have std == small, dead channel/noise
                ADV=0.0005;
                AScans(:,val < ADV) = 0;
                
                %detect bad artifact AScans
                ABV=0.6;
                AScans(:,val > ABV)=0;
                
                if (sum(val > ABV)+sum(val==ACV)+sum(val < ADV)) > 0 
                    disp(['AScans dropped: ' 'const: ' num2str(sum(val==ACV)) ', dead: ' num2str(sum(val<ADV)) ', bad: ' num2str(sum(val > ABV)) ', mean-std: ' num2str(mean(val))]);
                else
                    disp(['AScans mean std.: ' num2str(mean(val)) ]);
                end
                
            end
            
            if CE_SF>SF %to same samplefreq as AScans
                if mod(size(CE,1),2) CE(end+1,:)=median(CE(:)); end
                if any(CE_SF/SF == [2,4,6,8,10,20])
                  CE=interpft(CE,size(CE,1)*(CE_SF\SF));
                else
                  CE=interp1((0:length(CE)-1)',CE,(0:CE_SF/SF:length(CE)-1)','PCHIP',0);
                end
            end
           
            if length(CE)<size(AScans,1) %padding, again
                CE(size(AScans,1))=median(CE(:));                
            end
            CE=mean(CE,2); %if CE_measured
            CE=repmat(CE,[1 size(AScans,2)]); %padding for AScan amount
            
            %%upsample DATA if requested
            if SF<upsampling
              AScans=interpft(AScans,size(AScans,1)*ceil(upsampling/SF));
              CE=interpft(CE,size(CE,1)*ceil(upsampling/SF)); SF=upsampling; %update SF afterwards
            end
            
            %system time delay
            t0=(0:1/ SF:( size(AScans,1) -1)*1/ SF);
            t=t0+eoffset;
            if eoffset~=0
                if rem(eoffset,(1/SF)) == 0
                  AScans = circshift(AScans,fix(eoffset/(1/SF)));
                  %AScans(1:fix(eoffset/(1/SF)),:)=0;
                else
                  for i=1:size(AScans,2)
                    AScans(:,i)=interp1(t,AScans(:,i),t0,'PCHIP',mean(AScans(:,i)));
                  end
                end
            end
            %figure; plot(t0,data,t0,data_i);
           
            %%%recover amplitude
            AScans=AScans./repmat(Amplification',[size(AScans,1) 1]);
            
            if matchedFiltering==1
                %%%matched filtered bandwith restoring
                AScans=fft(AScans).*conj((fft(CE)+20*eps)./ abs((fft(CE)+20*eps)));
                AScans(1,:) = 0; %%DC free
                AScans=ifft(AScans);
            end
             
            if envelopeImaging==1
                AScans=abs(hilbert(AScans));
            end
                        
            for rT=1:numTAS
                for rE=1:numRec
                    %%angle filtering
                    inbetweenAngle3D=180*angleBetweenVectors(squeeze(geomRec(rE,:,rT))-middlePoint,squeeze(geomEmit(eE,:,eT))-middlePoint)/pi;
                    if(~isValidEmRecComb(angleEmRecComb,squeeze(geomRec(rE,:,rT))-middlePoint,squeeze(geomEmit(eE,:,eT))-middlePoint)), disp(['supressed: ' num2str([eT rT inbetweenAngle3D])]); continue; end
                    
                    i=i+1; %used ascans
                    soundvelocity=soundSpeed((TASTemperature(1,rT,Mp)+TASTemperature(1,eT,Mp))/2); %mean temp between
                                  
                    %selected Ascan
                    Data=AScans(:,find(receiverIndices==rE & TASIndices==rT));
                    %SF=10e6; %%from org Data -> will be changed
                    %figure(3); plot(Data)
                    if isempty(Data) continue; disp('warning: empty data'); end %%%early exit
                                       
                    if visualization>1
                        figure(1);
                        plot(Data); %visualization
                        title(sprintf('TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                    end                    
                
                    imagPos=imagstart;
                    if SAFT_MEX==1 || SAFT_MEX==2 %SAFT MEX
                        %if visualization>1 tic, end
                        imag=addsig2vol_3(circshift(Data,+0),single(imagPos'),single(squeeze(geomRec(rE,:,rT)))',single(squeeze(geomEmit(eE,:,eT)))',single(soundvelocity),single(imagRes),single(1/SF),uint32(imagXYZ'),imag);
                        %if visualization>1 t=toc, (t\numel(imagPosAll))/1024^2, end  % kVoxel/s
                    end
                    if SAFT_MEX==0 || SAFT_MEX==2 %%MATLAB SAFT
                        %%SAFT MATLAB NAIVE
                        %imag=SAFT_MATLAB(Data,imagstart,squeeze(geomRec(rE,:,rT)),squeeze(geomEmit(eE,:,eT)),soundvelocity,imagRes,SF,imagXYZ,imag);
                        
                        %%SAFT MATLAB_vec
                        %tic
                        %Data(length(Data)*3)=0; %padding to maximum possible size (for index access)
                        imag=SAFT_MATLAB_vec(Data,imagPosAll,squeeze(geomRec(rE,:,rT)),squeeze(geomEmit(eE,:,eT)),soundvelocity,imagRes,SF,imagXYZ,imag);
                        %t=toc, (t\numel(imagPosAll))/1024^2    % kVoxel/s
                    end
                    
                    if visualization>1
                        figure(2);
                        plot3(reshape(geomRec(:,1,:),[numTAS*numRec 1]),reshape(geomRec(:,2,:),[numTAS*numRec 1]),reshape(geomRec(:,3,:),[numTAS*numRec 1]),'.b',reshape(geomEmit(:,1,:),[numTAS*numEmit 1]),reshape(geomEmit(:,2,:),[numTAS*numEmit 1]),reshape(geomEmit(:,3,:),[numTAS*numEmit 1]),'.g',squeeze(geomRec(rE,1,rT)),squeeze(geomRec(rE,2,rT)),squeeze(geomRec(rE,3,rT)),'or',squeeze(geomEmit(eE,1,eT)),squeeze(geomEmit(eE,2,eT)),squeeze(geomEmit(eE,3,eT)),'or');
                        hold on
                        scatter3(meshx(:), meshy(:), meshz(:), 5.*ones([numel(meshy) 1]), reshape(imag(1:downsampling:end,1:downsampling:end,1:downsampling:end)',[numel(imag(1:downsampling:end,1:downsampling:end,1:downsampling:end)) 1]));
                        hold off
                        
                        set(gca,'ZDir','reverse'); title(sprintf('TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        drawnow;
                    end
                    
                    if visualization>0 && SAFT_MEX<2
                        %figure(3); imagesc(imagPosX,imagPosY,sqrt(abs(imag(:,:,1)))); title(sprintf('SQRT dedynamic SAFT: TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        figure(3); imagesc(imagPosX,imagPosY,((imag(:,:,1)))); title(sprintf('SAFT: TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        colorbar; drawnow;
                    else
                        figure(4); imagesc(imagPosX,imagPosY,imag2(:,:,1)); title(sprintf('SAFT MATLAB, TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        figure(3); imagesc(imagPosX,imagPosY,imag(:,:,1)); title(sprintf('SAFT MEX, TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        figure(5); imagesc(imagPosX,imagPosY,imag2(:,:,1)./max(max(imag2(:,:,1)))-imag(:,:,1)./max(max(imag(:,:,1)))); title(sprintf('DIFF SAFT, TAS-Receiver %i, Receiver ele. %i, TAS-Emitter %i, Emitter ele. %i, Movement position %i',rT, rE, eT, eE, Mp));
                        drawnow;
                    end
                    
                    %if SAFT_MEX~=1
                    %    imagsum=imagsum+imag;
                    %else
                    %     imagsum=imag;
                    %end
                end
                %profile report
                
            end
        end
        
        if image_save==1
            if visualization>0 & SAFT_MEX<2
                saveas(3,sprintf('.%cSAFT-MATLAB_20delay_ASCAN%06i_TR%i_RE%i_TE%i_EE%i_MP%i%s',filesep,i, rT, rE, eT, eE, Mp,'.png'),'png');
            else
                saveas(4,sprintf('.%cSAFT-MEX2_ASCAN%06i_TR%i_RE%i_TE%i_EE%i_MP%i%s',filesep,i, rT, rE, eT, eE, Mp,'.png'),'png');
                saveas(5,sprintf('.%cSAFT-DIFF2_ASCAN%06i_TR%i_RE%i_TE%i_EE%i_MP%i%s',filesep,i, rT, rE, eT, eE, Mp,'.png'),'png');
            end
        end
    end
end

%save('workspace.mat');
end

function imag=SAFT_MATLAB(Data,imagstart,geomRec,geomEmit,soundvelocity,imagRes,SF,imagXYZ,imag)
%%%naive reference implementaion
imagPos=imagstart;
for x=1:imagXYZ(1)
    
    imagPos(2)= imagstart(2);
    for y=1:imagXYZ(2)
        
        imagPos(3)= imagstart(3);
        for z=1:imagXYZ(3)
            
            dist=sqrt(sum((imagPos-geomRec).^2));
            dist=dist+sqrt(sum((imagPos-geomEmit).^2));
            %%+1 for matlab 1 indexing
            try imag(x,y,z)=imag(x,y,z)+Data(1+round(SF*(dist/soundvelocity))); catch, round(SF*(dist/soundvelocity)); end %%% out of array -> continue
            
            imagPos(3)= imagPos(3)+ imagRes;
        end
        imagPos(2)= imagPos(2)+ imagRes;
    end
    imagPos(1)= imagPos(1)+ imagRes;
end
end

function imag=SAFT_MATLAB_vec(Data,imagPosAll,geomRec,geomEmit,soundvelocity,imagRes,SF,imagXYZ,imag)
%vectorized MATLAB reference implementation, with Divide and conquer
try
    %Data(length(Data)*3)=0; %padding to maximum possible size (for index access)
    
    gR=repmat(shiftdim(geomRec,-2),imagXYZ);
    gE=repmat(shiftdim(geomEmit,-2),imagXYZ);
    
    dist=sqrt(sum((gR-imagPosAll).^2,4)) + sqrt(sum((gE-imagPosAll).^2,4));
    idx=1+round(SF*dist/soundvelocity); %+1 for matlab 1 indexing
    
    %%data map
    idx(idx>size(Data,1))=size(Data,1); %set to last value....FIXME
    if all(idx==size(Data,1)) disp('Warning: everything outside view window'); return; end
    imag=imag+Data(idx);
catch %recursive imagesize reduction on OOM
    clear dist idx gR gE
    %%separating along third dimension 
    imag1=SAFT_MATLAB_vec(Data,imagPosAll(:,:,1:floor(size(imagPosAll,3)/2),:),geomRec,geomEmit,soundvelocity,imagRes,SF,[imagXYZ(1:2) ,floor(size(imagPosAll,3)/2)],imag(:,:,1:floor(size(imagPosAll,3)/2)));
    imag2=SAFT_MATLAB_vec(Data,imagPosAll(:,:,1+floor(size(imagPosAll,3)/2):end,:),geomRec,geomEmit,soundvelocity,imagRes,SF,[imagXYZ(1:2) ,length(imag(1,1,1+floor(size(imagPosAll,3)/2):end))],imag(:,:,1+floor(size(imagPosAll,3)/2):end));
    
    imag=cat(3,imag1,imag2);
end

%catch %recursive
%Data(length(Data)*2)=0; %padding to maximum possible size (for index access)
%imag=SAFT_MATLAB_vec(Data,imagPosAll,geomRec,geomEmit,soundvelocity,imagRes,SF,imagXYZ,imag)
%end

end

function b = isValidEmRecComb(angleLimit,s,r)
b = (angleBetweenVectors(s,r)<=angleLimit./180*pi);
end

function a = angleBetweenVectors(v1,v2)
    %k=zeros(size(v1,1),size(v1,2));
    %for i=1:size(v1,1) for j=1:size(v1,2)
    %        k(i,j)=10.*eps+norm(squeeze(v1(i,j,:))).*norm(squeeze(v2(i,j,:)));
    %end,end
    v1n=v1./repmat(sqrt(sum(v1.^2,2)),[1 size(v1,2)]);
    v2n=v2./repmat(sqrt(sum(v2.^2,2)),[1 size(v2,2)]);
    a = acos(sum(v1n.*v2n,2));
    if any(isnan(v1n)|isinf(v2n)|isnan(v2n)|isinf(v1n)) a=angleBetweenVectors(v1+rand(1,1).*10.*eps(v1),v2+rand(1,1).*10.*eps(v2)); end %recursion
end