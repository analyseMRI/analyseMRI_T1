%% Copyright as per 01/09/2021 by F.I. Kandi
% v6 includes cerebral cortex
%% Reset
clear all
close all
clc

%% Load patient
%IMPORTANT! Make a folder for each patient separately in the source
%directory with as foldername the patient code. In this folder place the
%brain.nii and aseg.nii files, without adapting their names.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Welcome to analyseMRI version 6! \n');
fprintf('\n');
fprintf('MRI T1-weighted scan analysis software for analyzing the Hippocampus, Amygdala and Cerebral Cortex \n');
fprintf('Usable on T1-weighted 1 mm isotropic voxel resolution scans\n');
fprintf('This version includes a comparison to an average healthy control group \n');
fprintf('\n');
fprintf('Written and developed by F.I. Kandi (copyright 01-01-2022) \n');
fprintf('For questions: analysemri@gmail.com \n');
fprintf('\n');
fprintf('This software is free to use \n');
fprintf('Please cite the developer and software when publishing having used this software \n');
fprintf('\n');
prompt_p = 'Enter patient number YY+C/D+XXX+T1:\n';
patient = input(prompt_p, 's');
if isempty(patient)
    return
end
prompt_s = 'Enter source directory in double quotations with backward slash, end with backward slash: \n';
source = input(prompt_s);
source = strrep(source,"\","/");
if isempty(source)
    return
end %omzetten naar backward
prompt_t = 'Enter target directory in double quotations with backward slash, end with backward slash: \n';
target = input(prompt_t);
if isempty(target)
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

source = strrep(source,"\","/");
brain = niftiread(source+patient+"/brain.nii");
aseg = niftiread(source+patient+"/aseg.nii"); 

%create final result directory
mkdir(target+patient+"\Histo");
mkdir(target+patient+"\Box");

%Atlas values per structure 
LH = 17; %17  Left-Hippocampus 
RH = 53; %53  Right-Hippocampus   
LA = 18; %18  Left-Amygdala   
RA = 54; %54  Right-Amygdala
LCWM = 2; %2   Left-Cerebral-White-Matter 
RCWM = 41; %41  Right-Cerebral-White-Matter
LCX = 3; %3   Left-Cerebral-Cortex
RCX = 42; %42  Right-Cerebral-Cortex

%structures = {'LH','RH','LA','RA','LCWM','RCWM','LCX','RCX'};
structures = {17,53,18,54,2,41,3,42};

for count_str = 1:length(structures)
    STRUCTURE = structures{count_str};

    if STRUCTURE == LH
        STNAME = "Left-Hippocampus";
    elseif STRUCTURE == RH
          STNAME = "Right-Hippocampus";
    elseif STRUCTURE == LA
          STNAME = "Left-Amygdala";
    elseif STRUCTURE == RA
          STNAME = "Right-Amygdala";
    elseif STRUCTURE == LCWM
          STNAME = "Left-Cerebral-White-Matter";
    elseif STRUCTURE == RCWM
          STNAME = "Right-Cerebral-White-Matter";
    elseif STRUCTURE == LCX
          STNAME = "Left-Cerebral-Cortex";
    elseif STRUCTURE == RCX
          STNAME = "Right-Cerebral-Cortex";
    end


    %LOADING IMAGES, CREATING MASK, APPLYING MASK
    [m,n,o] = size(brain); %get size of 3D image 
    aseg_mask = uint8(zeros(m,n,o)); %create empty 3D mask 
    masked_im = uint8(zeros(m,n,o)); %create empty 3D masked image 

    for s = 1:1:m %create mask, pick out only hippocampus
        for t = 1:1:n
            for u = 1:1:o
                if aseg(s,t,u) == STRUCTURE %sort out voxels of specific brian structure
                    aseg_mask(s,t,u) = 1; %make them 1
                else
                    aseg_mask(s,t,u) = 0; %make rest of brain 0
                end
            end
        end
    end

    for x = 1:1:o %applying mask
        masked_im(:,:,x) = brain(:,:,x).*aseg_mask(:,:,x);
    end

    %EVALUATION OF WHOLE BRAIN FORLOOP
    cl_slice = uint8(zeros(m,n));
    intensity = [];
    slice_volumes = [];
    for x = 1:o %o is number of z-slices i.e. max z coordinate
        input = masked_im(:,:,x); %MUST BE uint8
        if nansum(input) == 0 
            continue
        elseif sum(input) == 0
            continue
        else
            [mn, md, vr, st, areapix, areamm, volume_slice, intensities] = slice_analysis(input, 1, 1);
            means(x) = mn;
            modes(x) = md;
            variances(x) = vr;
            stdeviations(x) = st;
            pixelareas(x) = areapix;
            mmareas(x) = areamm;
            slice_volumes(x) = volume_slice;
            intensity = [intensity; intensities]; %column matrix of all pixel intensities
        end
    end
    structure_volume(count_str) = sum(slice_volumes);
    avg_intensity(count_str) = mean(intensity); 
    mode_intensity(count_str) = mode(intensity);
    var_intensity(count_str) = var(intensity);
    stdev_intensity(count_str) = std(intensity);

    figure(5)
    histogram(intensity);
    title("Intensity histogram, " + STNAME + ", Volume = " + structure_volume(count_str) + " mm3, Mean = " + avg_intensity(count_str), 'FontSize', 10);
    xlabel('Intensity range 0-255'); %0 is black, 255 is white
    ylabel('Frequency');

    figure(6)
    boxplot(intensity);
    title("Intensity boxplot, " + STNAME + ", Volume = " + structure_volume(count_str) + " mm3, Mean = " + avg_intensity(count_str), 'FontSize', 10);
    ylabel('Intensity range 0-255');


    folderhisto = target+patient+"\Histo";
    folderbox = target+patient+"\Box";
    basehisto = sprintf("Histo_"+STNAME+"_"+patient);
    basebox = sprintf("Box_"+STNAME+"_"+patient);
    fullhisto = fullfile(folderhisto, basehisto);
    fullbox = fullfile(folderbox, basebox);
    saveas(figure(5), fullhisto, 'png')
    saveas(figure(5), fullhisto, 'eps')
    saveas(figure(6), fullbox, 'png')
    saveas(figure(6), fullbox, 'eps')

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% RELATIVE VOLUMES AND INTENSITIES
analytics_names(1) = "LH_volume";
analytics_names(2) = "LH_volume_rel";
analytics_names(3) = "LH_mean_intensity";
analytics_names(4) = "LH_mean_int_rel";
analytics_names(5) = "LH_mode_intensity";
analytics_names(6) = "LH_mode_int_rel";
analytics_names(7) = "LH_var_intensity"; 
analytics_names(8) = "LH_stdev_intensity";

analytics_names(9) = "RH_volume";
analytics_names(10) = "RH_volume_rel";
analytics_names(11) = "RH_mean_intensity";
analytics_names(12) = "RH_mean_int_rel";
analytics_names(13) = "RH_mode_intensity";
analytics_names(14) = "RH_mode_int_rel";
analytics_names(15) = "RH_var_intensity"; 
analytics_names(16) = "RH_stdev_intensity";

analytics_names(17) = "LA_volume";
analytics_names(18) = "LA_volume_rel";
analytics_names(19) = "LA_mean_intensity";
analytics_names(20) = "LA_mean_int_rel";
analytics_names(21) = "LA_mode_intensity";
analytics_names(22) = "LA_mode_int_rel";
analytics_names(23) = "LA_var_intensity"; 
analytics_names(24) = "LA_stdev_intensity";

analytics_names(25) = "RA_volume";
analytics_names(26) = "RA_volume_rel";
analytics_names(27) = "RA_mean_intensity";
analytics_names(28) = "RA_mean_int_rel";
analytics_names(29) = "RA_mode_intensity";
analytics_names(30) = "RA_mode_int_rel";
analytics_names(31) = "RA_var_intensity"; 
analytics_names(32) = "RA_stdev_intensity";

analytics_names(33) = "LCX_volume";
analytics_names(34) = "LCX_volume_rel";
analytics_names(35) = "LCX_mean_intensity";
analytics_names(36) = "LCX_mean_int_rel";
analytics_names(37) = "LCX_mode_intensity";
analytics_names(38) = "LCX_mode_int_rel";
analytics_names(39) = "LCX_var_intensity"; 
analytics_names(40) = "LCX_stdev_intensity";

analytics_names(41) = "RCX_volume";
analytics_names(42) = "RCX_volume_rel";
analytics_names(43) = "RCX_mean_intensity";
analytics_names(44) = "RCX_mean_int_rel";
analytics_names(45) = "RCX_mode_intensity";
analytics_names(46) = "RCX_mode_int_rel";
analytics_names(47) = "RCX_var_intensity"; 
analytics_names(48) = "RCX_stdev_intensity";


analytics = [];
analytics(1) = structure_volume(1);
analytics(2) = structure_volume(1)/structure_volume(5);
analytics(3) = avg_intensity(1);
analytics(4) = avg_intensity(1)/avg_intensity(5);
analytics(5) = mode_intensity(1);
analytics(6) = mode_intensity(1)/mode_intensity(5);
analytics(7) = var_intensity(1);
analytics(8) = stdev_intensity(1);

analytics(9) = structure_volume(2);
analytics(10) = structure_volume(2)/structure_volume(6);
analytics(11) = avg_intensity(2);
analytics(12) = avg_intensity(2)/avg_intensity(6);
analytics(13) = mode_intensity(2);
analytics(14) = mode_intensity(2)/mode_intensity(6);
analytics(15) = var_intensity(2);
analytics(16) = stdev_intensity(2);

analytics(17) = structure_volume(3);
analytics(18) = structure_volume(3)/structure_volume(5);
analytics(19) = avg_intensity(3);
analytics(20) = avg_intensity(3)/avg_intensity(5);
analytics(21) = mode_intensity(3);
analytics(22) = mode_intensity(3)/mode_intensity(5);
analytics(23) = var_intensity(3);
analytics(24) = stdev_intensity(3);

analytics(25) = structure_volume(4);
analytics(26) = structure_volume(4)/structure_volume(6);
analytics(27) = avg_intensity(4);
analytics(28) = avg_intensity(4)/avg_intensity(6);
analytics(29) = mode_intensity(4);
analytics(30) = mode_intensity(4)/mode_intensity(6);
analytics(31) = var_intensity(4);
analytics(32) = stdev_intensity(4);

analytics(33) = structure_volume(7);
analytics(34) = structure_volume(7)/structure_volume(5);
analytics(35) = avg_intensity(7);
analytics(36) = avg_intensity(7)/avg_intensity(5);
analytics(37) = mode_intensity(7);
analytics(38) = mode_intensity(7)/mode_intensity(5);
analytics(39) = var_intensity(7);
analytics(40) = stdev_intensity(7);

analytics(41) = structure_volume(8);
analytics(42) = structure_volume(8)/structure_volume(6);
analytics(43) = avg_intensity(8);
analytics(44) = avg_intensity(8)/avg_intensity(6);
analytics(45) = mode_intensity(8);
analytics(46) = mode_intensity(8)/mode_intensity(6);
analytics(47) = var_intensity(8);
analytics(48) = stdev_intensity(8);

analytics_relative = []; %contains only relative values and var/stdev, for internal use only
analytics_relative(1) = structure_volume(1)/structure_volume(5);
analytics_relative(2) = avg_intensity(1)/avg_intensity(5);
analytics_relative(3) = mode_intensity(1)/mode_intensity(5);
analytics_relative(4) = var_intensity(1);
analytics_relative(5) = stdev_intensity(1);

analytics_relative(6) = structure_volume(2)/structure_volume(6);
analytics_relative(7) = avg_intensity(2)/avg_intensity(6);
analytics_relative(8) = mode_intensity(2)/mode_intensity(6);
analytics_relative(9) = var_intensity(2);
analytics_relative(10) = stdev_intensity(2);

analytics_relative(11) = structure_volume(3)/structure_volume(5);
analytics_relative(12) = avg_intensity(3)/avg_intensity(5);
analytics_relative(13) = mode_intensity(3)/mode_intensity(5);
analytics_relative(14) = var_intensity(3);
analytics_relative(15) = stdev_intensity(3);

analytics_relative(16) = structure_volume(4)/structure_volume(6);
analytics_relative(17) = avg_intensity(4)/avg_intensity(6);
analytics_relative(18) = mode_intensity(4)/mode_intensity(6);
analytics_relative(19) = var_intensity(4);
analytics_relative(20) = stdev_intensity(4);

analytics_relative(21) = structure_volume(7)/structure_volume(5);
analytics_relative(22) = avg_intensity(7)/avg_intensity(5);
analytics_relative(23) = mode_intensity(7)/mode_intensity(5);
analytics_relative(24) = var_intensity(7);
analytics_relative(25) = stdev_intensity(7);

analytics_relative(26) = structure_volume(8)/structure_volume(6);
analytics_relative(27) = avg_intensity(8)/avg_intensity(6);
analytics_relative(28) = mode_intensity(8)/mode_intensity(6);
analytics_relative(29) = var_intensity(8);
analytics_relative(30) = stdev_intensity(8);

xlswrite(target+patient+"\analytics_"+patient+".xlsx", analytics)
xlswrite(target+patient+"\analytics_colnames.xlsx", analytics_names)

% CONTROL VECTOR
control_relative = double([0.01825761,0.592941526,0.512272727,117.0933763,10.7904738,0.01838217,0.609363245,0.530454545,109.0385437,10.43238184,0.007211185,0.610994006,0.547727273,79.14366076,8.865790691,0.008117502,0.631962634,0.564545455,81.14798359,8.966254197,1.049148219,0.586258847,0.513636364,119.680083,10.93325958,1.047135667,0.588255327,0.514545455,118.9627604,10.9005493]);

% display results
if analytics_relative(1) < control_relative(1)
    fprintf('Left Hippocampus has smaller volume than healthy control \n');
elseif analytics_relative(1) > control_relative(1)
    fprintf('Left Hippocampus has larger volume than healthy control \n');
end 
if analytics_relative(2) < control_relative(2)
    fprintf('Left Hippocampus has lower relative intensity mean than healthy control \n');
elseif analytics_relative(2) > control_relative(2)
    fprintf('Left Hippocampus has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(3) < control_relative(3)
    fprintf('Left Hippocampus has lower relative intensity mode than healthy control \n');
elseif analytics_relative(3) > control_relative(3)
    fprintf('Left Hippocampus has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(4) < control_relative(4)
    fprintf('Left Hippocampus has lower intensity variance than healthy control \n');
elseif analytics_relative(4) > control_relative(4)
    fprintf('Left Hippocampus has higher intensity variance than healthy control \n');
end 
if analytics_relative(5) < control_relative(5)
    fprintf('Left Hippocampus has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(5) > control_relative(5)
    fprintf('Left Hippocampus has higher intensity standard deviation than healthy control \n');
end 

if analytics_relative(6) < control_relative(6)
    fprintf('Right Hippocampus has smaller volume than healthy control \n');
elseif analytics_relative(6) > control_relative(6)
    fprintf('Right Hippocampus has larger volume than healthy control \n');
end 
if analytics_relative(7) < control_relative(7)
    fprintf('Right Hippocampus has lower relative intensity mean than healthy control \n');
elseif analytics_relative(7) > control_relative(7)
    fprintf('Right Hippocampus has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(8) < control_relative(8)
    fprintf('Right Hippocampus has lower relative intensity mode than healthy control \n');
elseif analytics_relative(8) > control_relative(8)
    fprintf('Right Hippocampus has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(9) < control_relative(9)
    fprintf('Right Hippocampus has lower intensity variance than healthy control \n');
elseif analytics_relative(9) > control_relative(9)
    fprintf('Right Hippocampus has higher intensity variance than healthy control \n');
end 
if analytics_relative(10) < control_relative(10)
    fprintf('Right Hippocampus has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(10) > control_relative(10)
    fprintf('Right Hippocampus has higher intensity standard deviation than healthy control \n');
end 
    
if analytics_relative(11) < control_relative(11)
    fprintf('Left Amygdala has smaller volume than healthy control \n');
elseif analytics_relative(11) > control_relative(11)
    fprintf('Left Amygdala has larger volume than healthy control \n');
end 
if analytics_relative(12) < control_relative(12)
    fprintf('Left Amygdala has lower relative intensity mean than healthy control \n');
elseif analytics_relative(12) > control_relative(12)
    fprintf('Left Amygdala has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(13) < control_relative(13)
    fprintf('Left Amygdala has lower relative intensity mode than healthy control \n');
elseif analytics_relative(13) > control_relative(13)
    fprintf('Left Amygdala has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(14) < control_relative(14)
    fprintf('Left Amygdala has lower intensity variance than healthy control \n');
elseif analytics_relative(14) > control_relative(14)
    fprintf('Left Amygdala has higher intensity variance than healthy control \n');
end 
if analytics_relative(15) < control_relative(15)
    fprintf('Left Amygdala has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(15) > control_relative(15)
    fprintf('Left Amygdala has higher intensity standard deviation than healthy control \n');
end 

if analytics_relative(16) < control_relative(16)
    fprintf('Right Amygdala has smaller volume than healthy control \n');
elseif analytics_relative(16) > control_relative(16)
    fprintf('Right Amygdala has larger volume than healthy control \n');
end 
if analytics_relative(17) < control_relative(17)
    fprintf('Right Amygdala has lower relative intensity mean than healthy control \n');
elseif analytics_relative(17) > control_relative(17)
    fprintf('Right Amygdala has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(18) < control_relative(18)
    fprintf('Right Amygdala has lower relative intensity mode than healthy control \n');
elseif analytics_relative(18) > control_relative(18)
    fprintf('Right Amygdala has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(19) < control_relative(19)
    fprintf('Right Amygdala has lower intensity variance than healthy control \n');
elseif analytics_relative(19) > control_relative(19)
    fprintf('Right Amygdala has higher intensity variance than healthy control \n');
end 
if analytics_relative(20) < control_relative(20)
    fprintf('Right Amygdala has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(20) > control_relative(20)
    fprintf('Right Amygdala has higher intensity standard deviation than healthy control \n');
end 

if analytics_relative(21) < control_relative(21)
    fprintf('Left Cerebral Cortex has smaller volume than healthy control \n');
elseif analytics_relative(21) > control_relative(21)
    fprintf('Left Cerebral Cortex has larger volume than healthy control \n');
end 
if analytics_relative(22) < control_relative(22)
    fprintf('Left Cerebral Cortex has lower relative intensity mean than healthy control \n');
elseif analytics_relative(22) > control_relative(22)
    fprintf('Left Cerebral Cortex has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(23) < control_relative(23)
    fprintf('Left Cerebral Cortex has lower relative intensity mode than healthy control \n');
elseif analytics_relative(23) > control_relative(23)
    fprintf('Left Cerebral Cortex has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(24) < control_relative(24)
    fprintf('Left Cerebral Cortex has lower intensity variance than healthy control \n');
elseif analytics_relative(24) > control_relative(24)
    fprintf('Left Cerebral Cortex has higher intensity variance than healthy control \n');
end 
if analytics_relative(25) < control_relative(25)
    fprintf('Left Cerebral Cortex has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(25) > control_relative(25)
    fprintf('Left Cerebral Cortex has higher intensity standard deviation than healthy control \n');
end 

if analytics_relative(26) < control_relative(26)
    fprintf('Right Cerebral Cortex has smaller volume than healthy control \n');
elseif analytics_relative(26) > control_relative(26)
    fprintf('Right Cerebral Cortex has larger volume than healthy control \n');
end 
if analytics_relative(27) < control_relative(27)
    fprintf('Right Cerebral Cortex has lower relative intensity mean than healthy control \n');
elseif analytics_relative(27) > control_relative(27)
    fprintf('Right Cerebral Cortex has higher relative intensity mean than healthy control \n');
end 
if analytics_relative(28) < control_relative(28)
    fprintf('Right Cerebral Cortex has lower relative intensity mode than healthy control \n');
elseif analytics_relative(28) > control_relative(28)
    fprintf('Right Cerebral Cortex has higher relative intensity mode than healthy control \n');
end 
if analytics_relative(29) < control_relative(29)
    fprintf('Right Cerebral Cortex has lower intensity variance than healthy control \n');
elseif analytics_relative(29) > control_relative(29)
    fprintf('Right Cerebral Cortex has higher intensity variance than healthy control \n');
end 
if analytics_relative(30) < control_relative(30)
    fprintf('Right Cerebral Cortex has lower intensity standard deviation than healthy control \n');
elseif analytics_relative(30) > control_relative(30)
    fprintf('Right Cerebral Cortex has higher intensity standard deviation than healthy control \n');
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright as per 01/09/2021 by F.I. Kandi

function [mn, md, vr, st, areapix, areamm, volume_slice, intensities] = slice_analysis(masked_im, pixelarea, slicethickness)

    % masked_image in uint8, area in mm2 and thickness in mm
    [a,b]=size(masked_im);
    masked_image = double(masked_im);
    mask = zeros(a,b);
    for x=1:1:a
        for y=1:1:b
            if masked_image(x,y)==0 %make masked out areas nan to not disturb correct grayscale analysis
                mask(x,y)=nan;
            else
                mask(x,y)=masked_image(x,y);
            end
        end
    end
 
    temp_denan = (mask(~isnan(mask))); %remove nan from double
    mn = mean(temp_denan); %mean
    md = mode(temp_denan); %mode (most frequent)
    vr = var(temp_denan); %variance (average of the squared differences from the mean)
    st = std(temp_denan); %standard deviation (square root of variance)
    fprintf('Intensity mean = %4.2f\n', mn);
    fprintf('Intensity mode = %4.2f\n', md);
    fprintf('Intensity variance = %4.2f\n', vr);
    fprintf('Intensity st. deviation = %4.2f\n', st);
    tfnan = isnan(mask); %find nan entries
    areapix = 0;
    for x=1:1:a
        for y=1:1:b
            if tfnan(x,y)==0 %find non nan entries
                areapix=areapix+1; %count non nan entries
            end
        end
    end
    areamm=areapix*pixelarea;
    volume_slice = areapix*pixelarea*slicethickness;
    fprintf('Area (pixels) = %4.2f\n', areapix);
    fprintf('Area (mm) = %4.2f\n', areamm);
    fprintf('Volume (mm2) = %4.2f\n', volume_slice);
    
    intensities=temp_denan; %column of all pixel intensities
end
