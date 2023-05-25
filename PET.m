
clc
clear
close all

%% Shahryar Ebrahimi
%% S.N. = 810196093
%% PET Single Subject Study
%% Loading ...

Volume  = zeros (87,65,26,12);

for i = 1:9  
    Volume(:,:,:,i) = analyze75read(['PET_motor\s8np01160em0' int2str(i) 'R']); 
end

for i = 10:12  
    Volume(:,:,:,i) = analyze75read(['PET_motor\s8np01160em' int2str(i) 'R']); 
end

Volume  = mat2gray(Volume);

%% Section 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% averaging

Volume_avg = zeros(87,65,26);

for i = 1:12

    Volume_avg = Volume_avg + Volume(:,:,:,i) ;
    
end

Volume_avg = Volume_avg/i ;

% Masking usign a suitable threshold

Mask  = zeros(87,65,26);
for i = 1:26  
    
    avg = Volume_avg(:,:,i);
    thr = graythresh(avg);
    
    msk = Mask(:,:,i);
    msk ( avg > thr )  = 1 ;
    Mask(:,:,i) = msk;

end

%  count number of non-zero pixels as intra-cerebral voxels

N   = nnz(Mask);
str = sprintf('\n* number of intra-cerebral voxels is %d', N);
fprintf(str);

% apply mask to the data

Volume_Mask = zeros (87,65,26,12);
for i  = 1:12
    Volume_Mask(:,:,:,i) = Mask .* Volume(:,:,:,i);
end

% plot intra-cerebral voxels Mask for odd slices (1-23)

figure('name','Intra-cerebral voxels','NumberTitle','off')

for i = 1:2:23
    subplot(3,4,(i+1)/2)
    imshow(Mask(:,:,i))
    title(['Slice #' int2str(i)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate average at rest

Rest_avg = zeros(87,65,26);

for i = 2:2:12  % at rest volumes
    Rest_avg = Rest_avg + Volume(:,:,:,i);
end

Rest_avg = Rest_avg/6;

% calculate average in task

Task_avg = zeros(87,65,26);

for i = 1:2:11  % at task volumes
    Task_avg = Task_avg + Volume(:,:,:,i);
end

Task_avg = Task_avg/6;

% calculate S_P

S_P = zeros(87,65,26);

for i=1:87
    for j=1:65
        for k=1:26
            
            Rest_norm = 0;
            Task_norm = 0;
            
            % calculate variance of pixel at rest
            
            for p = 2:2:12
                Rest_norm = Rest_norm + ( Volume(i,j,k,p) - Rest_avg(i,j,k) )^2;
            end
            
            % calculate variance of pixel in task
            
            for p = 1:2:11
                Task_norm = Task_norm + ( Volume(i,j,k,p) - Task_avg(i,j,k) )^2;
            end
            
            % calculate S_P
            
            S_P(i,j,k) = sqrt( ( Rest_norm + Task_norm )/(6+6-2) );
            
        end
    end
end

% calculate S

S  = sqrt( S_P.^2 * (1/6+1/6) );

% calculate t0 
t0 = (Task_avg - Rest_avg)./S;

% plot t0-value of voxels for odd slices(1-21) + 24th slice

figure('name','t0-value of voxels','NumberTitle','off')

for i = 1:2:21
    subplot(3,4,(i+1)/2)
    imshow(mat2gray(t0(:,:,i)).*Mask(:,:,i))
    title(['Slice #' int2str(i)])
end

% slice 24th
subplot(3,4,12)
imshow(mat2gray(t0(:,:,24)).*Mask(:,:,24))
title(['Slice #' int2str(24)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate threshold

p_value = 0.01;
thr     = abs( tinv(p_value,6+6-2) );

% display threshold value

str = sprintf('\n* Threshold amount of t_0 value is %d',thr);
fprintf(str);

% obtain avtive regions by thresholding

Active_map = zeros(87,65,26);
Active_map( t0 >= thr ) = 1;

% plot active voxels for odd slices(1-21) + 24th slice

figure('name','active voxels with p-value 0.01','NumberTitle','off')

for i = 1:2:21
    subplot(3,4,(i+1)/2)
    imshow(mat2gray(Active_map(:,:,i) + Volume_avg(:,:,i)).*Mask(:,:,i))
    title(['Slice #' int2str(i)])
end

% slice 24th
subplot(3,4,12)
imshow(mat2gray(Active_map(:,:,24) + Volume_avg(:,:,24)).*Mask(:,:,24))
title(['Slice #' int2str(24)])


%% Secton 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% will be available on Final_Report file in PDF format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 2  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% construct data matrix (X)

X  =  zeros(12,87*65*26 );

for i = 1:12
    X(i, :) = reshape( Volume(:,:,:,i), [1,87*65*26] );
end

%  calculate eigenvalues of X

X_eig  = eig (X * X');
X_eig  = sqrt(X_eig);
X_eig  = sort(X_eig,'descend');

% plot eigenvalues of data matrix

figure('name', 'Eigenvalues of data matrix','NumberTitle','off')
plot(X_eig,':*', 'color', 'b')
grid on
title('Eigenvalues of data matrix')

% normalize eigenvalues

X_eig = X_eig/sum(X_eig);

% plot normalized eigenvalues of data matrix

figure('name', 'Normalized eigenvalues of data matrix','NumberTitle','off')
plot(X_eig,':*', 'color', 'b')
grid on
title('Normalized eigenvalues of data matrix')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract eigenvariate & eigenimage related to first largest eigenvalue of X

[U,S,V] = svds(X,1);

% plot eigenvariate

figure('name', 'Eigenvariate related to the first largest eigenvalue of x matrix','NumberTitle','off')
i = 1:12;
plot(i,U,'-o' , 'color', 'b')
xlim([1 12])
grid on
title('Eigenvariate related to the first largest eigenvalue of x matrix')

% extract eigenimage related to first largest eigenvalue

Im_zy  = zeros(87, 65, 26);
Im_z   = vec2mat(V,87);
Im_z   = Im_z';

for i = 1:87
    hlp          = vec2mat(Im_z(i,:), 65);
    Im_zy(i,:,:) = hlp';
end

% plot eigenimage for odd slices(1-21) + 24th slice

figure('name','Eigenimage','NumberTitle','off')

for i = 1:2:21
    subplot(3,4,(i+1)/2)
    imshow(mat2gray(Im_zy(:,:,i)).*Mask(:,:,i))
    title(['Slice #' int2str(i)])
end

% slice 24th
subplot(3,4,12)
imshow(mat2gray(Im_zy(:,:,24)).*Mask(:,:,24))
title(['Slice #' int2str(24)])

% clear unnecessary values

clear avg hlp i j k p Im_z msk str 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Part 4  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% will be available on Final_Report file in PDF format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THE  END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
