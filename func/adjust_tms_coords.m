function adjusted_mni_coords = adjust_tms_coords(mni_coords)

    gmFile='INVERSE_of_MNI152_T1_2mm_brain.nii';
    
    gray_msk=niftiread(gmFile);
    
    ind=find(gray_msk);
    
    info = niftiinfo(gmFile);
    transform = info.Transform.T;
    transform = transform';
    
    rob_inv = inv(transform);
    
    adjusted_mni_coords = zeros(size(mni_coords));
    
    for c = 1:length(mni_coords)
        x1=mni_coords(c,1);
        y1=mni_coords(c,2);
        z1=mni_coords(c,3);
        u1=[x1,y1,z1,1]';
        u2=rob_inv*u1;
        d=zeros(length(ind),1);
        for i=1:length(ind)
            [xx,yy,zz]=ind2sub(size(gray_msk),ind(i));
            d(i)=sqrt( (xx-u2(1))^2+(yy-u2(2))^2+(zz-u2(3))^2 );
        end
        [~,ind_min]=min(d);
        [xx,yy,zz]=ind2sub(size(gray_msk),ind(ind_min));
        Voxmm=2;
        
        data.xyz=[xx,yy,zz]*Voxmm;
        tmp_xyz = [data.xyz(1,:)];
        VoxDims=2; Cntr=[46,64,37]; NegateX=1;   %correct for MNI 2mm
        clearvars x y ;
        
        x=tmp_xyz/VoxDims;     %convert to voxel index
        y=x-Cntr;
        if NegateX
            y(1)=y(1)*-1;
        end
        y=y*VoxDims;
        
        MNI_out = round(y);
        adjusted_mni_coords(c,1) = MNI_out(1);
        adjusted_mni_coords(c,2) = MNI_out(2);
        adjusted_mni_coords(c,3) = MNI_out(3);

    end
    
end