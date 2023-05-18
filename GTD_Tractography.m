function GTD_Tractography()
%Generate the circle data
posit_ic = 1;
R=0;
line_circles_predefine = zeros(401,2,5);
mask_location = zeros(100,100);
roi_position = zeros(1,3);
directions = zeros(1,3);
 for  curr_width = 1:16
         r=20+curr_width*0.55;
        theta=0:pi/200:2*pi;
        x=r*cos(theta)+30; 
        y=r*sin(theta)+30;    
        coor_s = round(curr_width);
       line_circles_predefine(:,1,coor_s) = x;
       line_circles_predefine(:,2,coor_s) = y;
%     plot(x,y,'-');
%       hold on;
 end
all_mat_directions = zeros(100,100,2);
start_points_roi = zeros(1,3);
end_points_roi = zeros(1,3);
   for  point_numb = 1:392    
       for  follor = 1:5
         for  line_numb = 1:16 
         curr_point_val = line_circles_predefine(point_numb,:,line_numb);
         curr_coor_x = round(curr_point_val(1,1));
         curr_coor_y = round(curr_point_val(1,2));
         if  curr_coor_x<=0
             curr_coor_x=1;
         end
         if  curr_coor_y<= 0
             curr_coor_y=1;
         end
          cur_dirction1 = line_circles_predefine(point_numb,:,line_numb);
          cur_dirction2=line_circles_predefine(point_numb+1,:,line_numb);   
          cur_dirction = cur_dirction1-cur_dirction2;
          leg=sqrt(sum(cur_dirction.^2));
          cur_dirction =cur_dirction/leg;
       if  posit_ic==1
         all_mat_directions(curr_coor_x,curr_coor_y,:) = -cur_dirction;   
       end
       if  posit_ic==0
             all_mat_directions(curr_coor_x,curr_coor_y,:) = cur_dirction;   
       end
         roi_position = [roi_position;[curr_coor_x,curr_coor_y,follor]];
         directions = [directions; cur_dirction,0];
         mask_location(curr_coor_x,curr_coor_y) =1;
       if   (point_numb< 50 )&&(point_numb> 1)        
           start_points_roi = [start_points_roi;[curr_coor_x,curr_coor_y],follor];
       end
    if   (point_numb > 340 )&&(point_numb <=392)        
           end_points_roi = [end_points_roi;[curr_coor_x,curr_coor_y],follor];
    end
         end
       end
   end
  if  posit_ic==1
     ror_seeds = start_points_roi;
  end
    if  posit_ic==0
        ror_seeds = end_points_roi;
    end
 [VectorPlotX,VectorPlotY]=meshgrid(1:100,1:100);
%quiver(VectorPlotX,VectorPlotY,all_mat_directions(:,:,1)',all_mat_directions(:,:,2)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=100; % the region
r0=30;
r1=20;% the inner radio of the circle
r2=30;% the outer radio of the circle
k=10;%The number of the circle
%The noise level as being R%
Slice=28;
l_1=0.002;
l_2=0.0006;
l_3=0.0006;
lamda=l_2*ones(I,I,3);
S0=35*ones(I,I,Slice);
b=[0,800*ones(1,6)];
G=[0 0 0;1 0 1;-1 0 1;0 1 1;0 1 -1; 1 1 0;-1 1 0];%G=[0,0;0,1;1,0;1,1];The gradient deirection
for i=2:7
    G(i,:)=G(i,:)./norm(G(i,:));
end
S=15*ones(I,I,Slice,7);
S_n=S;
FA=zeros(I,I,Slice);
VectorF=rand*ones(I,I,Slice,3);

for i=1:100
    for j=1:100
       if ((i-r0)^2+(j-r0)^2)<r2^2 &((i-r0)^2+(j-r0)^2)> r1^2;   
            alfa=atan((j-r0)/(i-r0));
            if i>=r0
                alfa=alfa;
            else
                alfa=pi+alfa;
            end
            v1=[-sin(alfa) cos(alfa) 0];
            lamda(i,j,1)=l_1;
            lamda(i,j,2)=l_2;
            alfa(i,j)=lamda(i,j,2);
            beta(i,j)=lamda(i,j,1)-alfa(i,j);
       else
           v1=[rand rand 0];
           v1=v1./norm(v1);
       end
       VectorF(i,j,13,:)=v1;
       ADCv=(lamda(i,j,1)+lamda(i,j,2)+lamda(i,j,2))/3;
       FAv=sqrt(1.5)*( sqrt((lamda(i,j,1)-ADCv).^2+(lamda(i,j,2)-ADCv).^2+(lamda(i,j,2)-ADCv).^2)./sqrt(lamda(i,j,1).^2+lamda(i,j,2).^2+lamda(i,j,2).^2) );
       FA(i,j,13)=FAv;
       v2(1)=-v1(2);
       v2(2)=v1(1);
       v2(3)=0;
       v3=[0 0 1];
       v3=v3./norm(v3);
       v=[ v1' v2' v3'];
       T=[lamda(i,j,1) 0 0;0 lamda(i,j,2) 0;0 0 lamda(i,j,3)];
       D=v*T*v';
       for k=1:size(G,1)
           logS(k)=log(S0(i,j,k))-b(k)*G(k,:)*D*G(k,:)';
       end
%      logS_n=logS+randn(size(logS))*(max(logS(2:7))-min(logS(2:7)))*R/100;
%      d=(max(logS(2:7))-min(logS(2:7)))*R/100
       S(i,j,13,:)=exp(logS);%compute the SiS
       %S_n(i,j,13,:)=exp(logS_n);         
    end
end
S_n=S+randn(size(S))*(R+20)/10;
S_n(S_n<0)=0.00001;
for i=1:7
    DTIdata(i).VoxelData =S_n(:,:,:,i); 
    DTIdata(i).Gradient = G(i,:);
    DTIdata(i).Bvalue=800;
end
DTIdata(1).Bvalue=0;
% Constants DTI
parametersDTI=[];
parametersDTI.BackgroundTreshold=0.05;
parametersDTI.WhiteMatterExtractionThreshold=0.10;
parametersDTI.textdisplay=true;
% Perform DTI calculation
[ADC,FA1,VectorF1,DifT]=DTI(DTIdata,parametersDTI);
%save('C:/Users/Administrator/Desktop/SDE-2021-8-26-paper_TMI/directions/DTI_directions.mat','S_n');
[VectorPlotX,VectorPlotY]=meshgrid(1:size(FA1,1),1:size(FA1,2));
quiver(VectorPlotX,VectorPlotY,VectorF1(:,:,13,1)',VectorF1(:,:,13,2)');
DTI_direc = VectorF1(:,:,13,:);
DTI_direc = squeeze(DTI_direc);
new_DTI_direc_correct = zeros(100,100,3);
numb=0;
ROI_Dirs  = zeros(0,3); 
ROI_positions = zeros(0,3);
for  curr_x_id = 1:100
    for  curr_y_id = 1:100
        cur_direc = squeeze(all_mat_directions(curr_x_id,curr_y_id,:));
        cur_direc_xx = cur_direc(1,1);
        cur_direc_yy = cur_direc(2,1);
        if  (cur_direc_xx~=0)||(cur_direc_yy~=0)
            curr_DTI_direc = squeeze(DTI_direc(curr_x_id,curr_y_id,:));
            curr_DTI_direc = curr_DTI_direc';
            curr_DTI_direc2 = curr_DTI_direc(1:2);
            leg=sqrt(sum(curr_DTI_direc2.^2));
            curr_DTI_direc2 =curr_DTI_direc2/leg;         
            multp_mat_val = cur_direc'*curr_DTI_direc2';
            if  multp_mat_val <0
                 curr_DTI_direc2 = -curr_DTI_direc2;
                 corner =  (acos(dot(cur_direc,curr_DTI_direc2)/(norm(cur_direc)*norm(curr_DTI_direc2))))*180/pi;                 
                 if  corner > 80
                 %    break;
                 end
                 new_DTI_direc_correct(curr_x_id,curr_y_id,:) = [curr_DTI_direc2,0];   
                 cur_roi_posi = [curr_x_id,curr_y_id,1];
                 ROI_positions = [ROI_positions;cur_roi_posi];                 
                 ROI_Dirs = [ROI_Dirs;curr_DTI_direc2,0];
            end
           if  multp_mat_val >=0
                corner =  (acos(dot(cur_direc,curr_DTI_direc2)/(norm(cur_direc)*norm(curr_DTI_direc2))))*180/pi;                 
                 if  corner > 80
                  %   break;
                 end
                 new_DTI_direc_correct(curr_x_id,curr_y_id,:) = [curr_DTI_direc2,0] ;   
                 cur_roi_posi = [curr_x_id,curr_y_id,1];
                 ROI_Dirs = [ROI_Dirs;curr_DTI_direc2,0];
                 ROI_positions = [ROI_positions;cur_roi_posi];
           end
        end
    end
end
quiver(VectorPlotX,VectorPlotY,VectorF1(:,:,13,1)',VectorF1(:,:,13,2)');
%quiver(VectorPlotX,VectorPlotY,all_mat_directions(:,:,1)',all_mat_directions(:,:,2)');
%quiver(VectorPlotX,VectorPlotY,new_DTI_direc_correct(:,:,1)',new_DTI_direc_correct(:,:,2)');
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROI_positions2 = ROI_positions;
ROI_positions2(:,3) = 2;
ROI_positions3 = ROI_positions;
ROI_positions3(:,3) = 3;
ROI_positions4 = ROI_positions;
ROI_positions4(:,3) = 4;
ROI_positions5 = ROI_positions;
ROI_positions5(:,3) = 5;
ROI_positions6 = ROI_positions;
ROI_positions6(:,3) = 6;
%%%%%%%%%%%
ROI_positions = [ROI_positions;ROI_positions2;ROI_positions3;ROI_positions4;ROI_positions5;ROI_positions6]
ROI_Dirs = [ROI_Dirs;ROI_Dirs;ROI_Dirs]
ROI_Dirs =[ROI_Dirs;ROI_Dirs];
ROIpositions = ROI_positions;
DirsROI = double(ROI_Dirs);
WeightedDirsROI = double(ROI_Dirs);
SDE_Order = 3;
step_size = 0.5;
str=['Preparing   fitting ...........'];
disp(str);
if  SDE_Order==3
   [A, fvalIntra1] = GetATernaryCubic(ROIpositions, DirsROI, WeightedDirsROI);
end
if  SDE_Order==4
   [A, fvalIntra1] = GetATernaryFifth(ROIpositions, DirsROI, WeightedDirsROI);
end
if  SDE_Order==5
   [A, fvalIntra3] = GetATernarySixth(ROIpositions, DirsROI, WeightedDirsROI);
end
if  SDE_Order==6
    [A, fvalIntra4] = GetATernarySeventh(ROIpositions, DirsROI, WeightedDirsROI);
end
if  SDE_Order==7
     [A, fvalIntra6] = GetATernaryEighth(ROIpositions, DirsROI, WeightedDirsROI);
end
Aall = cell(1,9); 
current_order = A;
Aall{1} = current_order; 
Aall{2} = current_order; 
Aall{3} = current_order; 
Aall{4} = current_order; 
Aall{5} = current_order; 
Aall{6} = current_order; 
Aall{7} = current_order; 
Aall{8} = current_order; 
Aall{9} = A; 

seeds=get_seed_wjq(ror_seeds);   
    Tracts = cell(1, 20000);
   for FiberNum = 1:size(seeds,1)
       disp(FiberNum);
        start_point = seeds(FiberNum,:);      
        if ~isempty(current_order) && any(current_order(:)~=0)
            FTrace = tracking_circle(Aall,start_point,step_size,mask_location);      
        end
        if ~isempty(FTrace)
            Tract = FTrace;       
        else
            Tract = [];
        end      
%%%%%%%%%%%%%%        
        intTract = round(Tract);     
        THRESH_VC = 126/step_size;
        if  size(Tract)<  THRESH_VC
        continue;
        end
        Tracts{FiberNum  } = Tract; 
   end
 %%%%%%%%%%%%%%%%%%%%%%%
 file_dir_all = 'C:/Users/Administrator/Desktop/SDE-2021-8-26-paper_TMI/CODE/github';
 peak_path = [file_dir_all,'/DTI_directions_peaks.mat'];
 cir_peaks = load(peak_path);
 mask_path = [file_dir_all,'/dti.nii.gz']; 
DWIdata = niftiread(mask_path); 
DWIinfo = niftiinfo(mask_path);
DWIaffine = (DWIinfo.Transform.T)'; 
MAT_NUIT = [1,0,0,-1;0,1,0,-0.7;0,0,1,-1;0,0,0,0];
Tracts = transform(Tracts,MAT_NUIT);
streamlines_saves = transform(Tracts, (DWIaffine));
save_img_vtk_copy.data = streamlines_saves;
save_img_vtk_copy.count = size(streamlines_saves,1);
save_img_vtk_copy.total_count = size(streamlines_saves,1);  
if  posit_ic==1
    path_saved = [file_dir_all,'/fiber/1Circle_fiber.tck'];
end
if  posit_ic==0
    path_saved = [file_dir_all,'/fiber/0Circle_fiber.tck'];
end
write_mrtrix_tracks(save_img_vtk_copy,path_saved)    ;
 
 
 
 
 
 
 
 
 
 