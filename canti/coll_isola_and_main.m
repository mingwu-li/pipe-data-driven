function [main_collsol,isola_collsol] = coll_isola_and_main(M,C,K,fnl,fFull,lambda,fext,phiend,os,newcoll)

if newcoll
    figure;
    DS = DynamicalSystem();
    set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
    set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
    set(DS.Options,'BaseExcitation',true);
    nmodes = size(M,1);
    DS.Omega = 0.95*abs(imag(lambda(3)));
    nCycles = 10;
    outdof = 1; % first modal coordinate
    kappas = [-1; 1];
    freqRange = [0.95 1.15]*abs(imag(lambda(3)));
    for k=1:numel(fFull)
        coeffs = [fext fext]/2;
        DS.add_forcing(coeffs, kappas, fFull(k));
        coll = cocoWrapper(DS, nCycles, outdof);
        set(coll.Options, 'dir_name', ['main_po_force_',num2str(k)],...
            'h_max', 2, 'NAdapt', 5, 'PtMX', 200, 'NSV',5);        
        coll.initialGuess = 'linear';
        coll.extract_FRC(freqRange);
    end
    % post-processing: from modal coordinates to physical coordinates
    main_collsol = cell(numel(fFull),1);
    for k=1:numel(fFull)
        runid = ['main_po_force_',num2str(k),'.FRC'];
        bd = coco_bd_read(runid);
        EPlab = coco_bd_labs(bd,'EP');
        EPlab = max(EPlab);
        omegas = zeros(EPlab,1);
        yends  = zeros(EPlab,1);
        for j=1:EPlab
            solj = po_read_solution(runid,j);
            yend = solj.xbp(:,1:nmodes)*phiend;
            omegas(j) = solj.p(1);
            yends(j)  = norm(yend,'inf');
        end
        sol = struct();
        sol.omega = omegas; sol.yend = yends;
        main_collsol{k} = sol;
    end
    save('dataVKcoll_post_main.mat','main_collsol')
else
    load('dataVKcoll_post_main.mat', 'main_collsol')
end
%% Isola part via collocation method
if newcoll
    figure;
    DS = DynamicalSystem();
    set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
    set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
    set(DS.Options,'BaseExcitation',true);
    DS.Omega = os;
    nCycles = 100;
    outdof = 1; % first modal coordinate
    kappas = [-1; 1];
    freqRange = [0.95 1.15]*abs(imag(lambda(3)));
    for k=1:numel(fFull)
        coeffs = [fext fext]/2;
        DS.add_forcing(coeffs, kappas, fFull(k));
        coll = cocoWrapper(DS, nCycles, outdof);
        set(coll.Options, 'dir_name', ['isola_po_force_',num2str(k)],...
            'NTST',20,'h_max', 10, 'NAdapt', 2, 'PtMX', 100, 'NSV',5,'bi_direct',false);        
        coll.initialGuess = 'forward';
        coll.extract_FRC(freqRange);
    end
    % post-processing: from modal coordinates to physical coordinates
    isola_collsol = cell(numel(fFull),1);
    for k=1:numel(fFull)
        runid = ['isola_po_force_',num2str(k),'.FRC'];
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
        isola_collsol{k} = sol; 
    end
    save('dataVKcoll_post_isola.mat','isola_collsol')
else
    load('dataVKcoll_post_isola.mat','isola_collsol')
end

end