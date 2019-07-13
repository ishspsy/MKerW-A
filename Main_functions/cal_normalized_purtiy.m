
n=660;
puri=[]; nnmi=[];  rrand=[];
for irep=1:1000
pur=[]; nnm=[]; rran=[];
for CC=2:30;
tt1 = randi([1 CC],1,n);  tt2 = randi([1 CC],1,n);
pur=[pur,purity(CC, tt1, tt2)];  nnm=[nnm,Cal_NMI(tt1, tt2)]; rran=[rran, RandIndex(tt1, tt2)];
end
puri=[puri;pur]; nnmi=[nnmi;nnm]; rrand=[rrand;rran];
end

norm_cls_purity_const=[mean(puri); std(puri)];
norm_cls_nmi_const=[mean(nnmi); std(nnmi)];
norm_cls_rand_const=[mean(rrand); std(rrand)];


save('norm_cls_number_const.mat', 'norm_cls_purity_const','norm_cls_nmi_const','norm_cls_rand_const')