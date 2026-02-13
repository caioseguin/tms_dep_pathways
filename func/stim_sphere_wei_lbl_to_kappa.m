function kappa = stim_sphere_wei_lbl_to_kappa(stim_sphere_wei_lbl)

    switch stim_sphere_wei_lbl
        case 'uni'
            kappa = 0;   
        case 'em0d05'
            kappa = -0.05;
        case 'em0d25'
            kappa = -0.25;
        case 'em0d5'
            kappa = -0.5;
        case 'em1'
            kappa = -1;
        case 'em1d5'
            kappa = -1.5;
        case 'em2'
            kappa = -2;
        otherwise
            error('Undefined stim_sphere_wei_lbl');
    end

end