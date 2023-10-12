function [bd1,bd2] = coll_SN_cont(M,C,K,fnl,fFull,omegaSpan,epsfSpan,fext,os,sncoll)

%% collocation of SN points
if sncoll
    DS = DynamicalSystem();
    set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
    set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex');
    set(DS.Options,'BaseExcitation',true);
    DS.Omega = os;
    nCycles = 10;
    outdof = 1; % first modal coordinate
    kappas = [-1; 1];
    coeffs = [fext fext]/2;
    DS.add_forcing(coeffs, kappas, fFull(3));
    coll  = cocoWrapper(DS, nCycles, outdof);  
    set(coll.Options, 'dir_name', 'postSN1', 'NSV', 10,'bi_direct',true);
    runid = 'isola_po_force_3.FRC'; bdtmp = coco_bd_read(runid);
    SNlab = coco_bd_labs(bdtmp,'SN');
    bd1 = coll.extract_SNFRC(runid, SNlab(1), omegaSpan, epsfSpan);
    set(coll.Options, 'dir_name', 'postSN2', 'NSV', 10);
    bd2 = coll.extract_SNFRC(runid, SNlab(2), omegaSpan, epsfSpan);
else
    bd1 = coco_bd_read('postSN1.SNFRC');
    bd2 = coco_bd_read('postSN2.SNFRC');
end

end