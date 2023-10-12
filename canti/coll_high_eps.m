function hig_collsol = coll_high_eps(M,C,K,fnl,fext,soleps,lambda,phiend,newcoll)


%% collocation validation
if newcoll
    figure;
    fFullsn = soleps.p(1);
    DS = DynamicalSystem();
    set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
    set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
    set(DS.Options,'BaseExcitation',true);
    nmodes = size(M,1);
    DS.Omega = 1.15*abs(imag(lambda(3)));
    nCycles = 10;
    outdof = 1; % first modal coordinate
    kappas = [-1; 1];
    freqRange = [0.95 1.15]*abs(imag(lambda(3)));
    coeffs = [fext fext]/2;
    DS.add_forcing(coeffs, kappas, fFullsn);
    coll = cocoWrapper(DS, nCycles, outdof);
    set(coll.Options, 'dir_name', 'high_po_force',...
        'h_max', 2, 'NAdapt', 5, 'PtMX', 200, 'NSV',5, 'ItMX', 15);        
    coll.initialGuess = 'linear';
    coll.extract_FRC(freqRange);
    % post-processing: from modal coordinates to physical coordinates
    runid = 'high_po_force.FRC';
    bd = coco_bd_read(runid);
    EPlab = coco_bd_labs(bd,'EP');
    EPlab = max(EPlab);
    omegas = zeros(EPlab,1);
    yends  = zeros(EPlab,1);
    stabs  = false(EPlab,1);
    for j=1:EPlab
        solj = po_read_solution(runid,j);
        yend = solj.xbp(:,1:nmodes)*phiend;
        omegas(j) = solj.p(1);
        yends(j)  = norm(yend,'inf');
        stabs(j)  = all(abs(solj.po_test.la)<1);
    end
    sol = struct();
    sol.omega = omegas; sol.yend = yends; sol.stabs = stabs;
    hig_collsol = sol;
    save('dataVKcoll_post_hig.mat','hig_collsol')
else
    load('dataVKcoll_post_hig.mat', 'hig_collsol')
end

end