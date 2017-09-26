nGrids=8;

% assumedSoln='constant';
% [error_phi0_n_FS_constant, order_phi_nMinus1_FS_constant]=converger_FS(assumedSoln,nGrids);
% [error_phi0_n_LS_constant, order_phi_nMinus1_LS_constant]=converger_LS(assumedSoln,nGrids);
% 
% assumedSoln='linear';
% [error_phi0_n_FS_linear, order_phi_nMinus1_FS_linear]=converger_FS(assumedSoln,nGrids);
% [error_phi0_n_LS_linear, order_phi_nMinus1_LS_linear]=converger_LS(assumedSoln,nGrids);
% 
% assumedSoln='quadratic';
% [error_phi0_n_FS_quadratic, order_phi_nMinus1_FS_quadratic]=converger_FS(assumedSoln,nGrids);
% [error_phi0_n_LS_quadratic, order_phi_nMinus1_LS_quadratic]=converger_LS(assumedSoln,nGrids);
% 
% assumedSoln='plus1Sqrt';
% [error_phi0_n_FS_plus1Sqrt, order_phi_nMinus1_FS_plus1Sqrt]=converger_FS(assumedSoln,nGrids);
% [error_phi0_n_LS_plus1Sqrt, order_phi_nMinus1_LS_plus1Sqrt]=converger_LS(assumedSoln,nGrids);

% reference solution
%%
% =================
% assumedSoln: constant
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_constant =   1.0e-13 *[...
   0.694513144194782
   0.761687275872297
   0.782914840154826
   0.788965484740937
   0.790013848361801
   0.790634961693368
   0.790438718455841
   0.788621155924496
];

order_phi_nMinus1_FS_constant =[...
  -0.133196796841732
  -0.039656594641704
  -0.011106797273031
  -0.001915755561086
  -0.001133809311344
   0.000358135298305
   0.003321203522507
];

% =================
% assumedSoln: constant
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_constant =   1.0e-13 *[...
   0.112228117813608
   0.112827285436508
   0.037772065111694
   0.039343271438244
   0.039291936301230
   0.011378431991337
   0.011370643298715
   0.008077778563184
];

order_phi_nMinus1_LS_constant =[...
  -0.007681825790841
   1.578724434669941
  -0.058797261970952
   0.001883658993058
   1.787931505097544
   0.000987882673690
   0.493283374614361
];
%%
% =================
% assumedSoln: linear
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_linear =[...
   0.048861548575326
   0.014391604851914
   0.003842496606856
   0.000980105451355
   0.000246338131562
   0.000061668161856
   0.000015422293693
   0.000003855902141
];

order_phi_nMinus1_FS_linear =[...
   1.763472107459269
   1.905111589342248
   1.971035101020718
   1.992297014744958
   1.998042215418898
   1.999508497337394
   1.999877004407445
];

% =================
% assumedSoln: linear
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_linear =   1.0e-13 *[...
   0.074257240236188
   0.047249088475604
   0.143580641605176
   0.125843669651456
   0.135971729298907
   0.156788099613461
   0.169156241165121
   0.452642638142466
];

order_phi_nMinus1_LS_linear =[...
   0.652245201271777
  -1.603502847163576
   0.190228603518118
  -0.111674077199480
  -0.205509338252473
  -0.109540346096509
  -1.420016082758704
];

%%
% =================
% assumedSoln: quadratic
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_quadratic =[...
   0.528435410891271
   0.155155848101193
   0.041512603380145
   0.010600167985160
   0.002665104766903
   0.000667238453540
   0.000166870090430
   0.000041721307520
];

order_phi_nMinus1_FS_quadratic =[...
   1.768009070358906
   1.902096760011719
   1.969462282578579
   1.991822975503822
   1.997917907248369
   1.999477043810538
   1.999869114288267
];

% =================
% assumedSoln: quadratic
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_quadratic =[...
   0.003582325361758
   0.000323354812612
   0.000023611148136
   0.000001548683201
   0.000000098048979
   0.000000006148216
   0.000000000384590
   0.000000000024048
];

order_phi_nMinus1_LS_quadratic =[...
   3.469706387034262
   3.775577980598491
   3.930354235224089
   3.981395632358030
   3.995262776263439
   3.998775409529993
   3.999325923258120
];

%%
% =================
% assumedSoln: plus1Sqrt
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_plus1Sqrt =[...
   0.007661883418932
   0.002481560333717
   0.000696777858899
   0.000180787262340
   0.000045651918675
   0.000011442186213
   0.000002862384315
   0.000000715711168
];

order_phi_nMinus1_FS_plus1Sqrt =[...
   1.626451543246898
   1.832476844935340
   1.946405747076788
   1.985545633219298
   1.996312767118894
   1.999073436593434
   1.999767989571842
];

% =================
% assumedSoln: plus1Sqrt
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_plus1Sqrt =   1.0e-03 *[...
   0.231940208089442
   0.024814379617639
   0.001929854481586
   0.000128936200543
   0.000008202699478
   0.000000514984548
   0.000000032222751
   0.000000002014445
];

order_phi_nMinus1_LS_plus1Sqrt =[...
   3.224504649159258
   3.684612414400402
   3.903762785470991
   3.974414792659945
   3.993497723792533
   3.998377548780968
   3.999625595867425
   ];
 
refinementRatio=2;
Tau=10;
for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
  gridMeshSize_n(iGrid)=Tau/J;
end

nSoln=4;
close all;
for iSoln=1:nSoln
  switch(iSoln)
    case 1 % constant
      soln='constant';
      error_phi0_n_FS=error_phi0_n_FS_constant;
      order_phi_nMinus1_FS=order_phi_nMinus1_FS_constant;
      error_phi0_n_LS=error_phi0_n_LS_constant;
      order_phi_nMinus1_LS=order_phi_nMinus1_LS_constant;
    case 2 % linear
      soln='linear';
      error_phi0_n_FS=error_phi0_n_FS_linear;
      order_phi_nMinus1_FS=order_phi_nMinus1_FS_linear;
      error_phi0_n_LS=error_phi0_n_LS_linear;
      order_phi_nMinus1_LS=order_phi_nMinus1_LS_linear;
    case 3 % quadratic
      soln='quadratic';
      error_phi0_n_FS=error_phi0_n_FS_quadratic;
      order_phi_nMinus1_FS=order_phi_nMinus1_FS_quadratic;
      error_phi0_n_LS=error_phi0_n_LS_quadratic;
      order_phi_nMinus1_LS=order_phi_nMinus1_LS_quadratic;
    case 4 % plus1Sqrt
      soln='plus1Sqrt';
      error_phi0_n_FS=error_phi0_n_FS_plus1Sqrt;
      order_phi_nMinus1_FS=order_phi_nMinus1_FS_plus1Sqrt;
      error_phi0_n_LS=error_phi0_n_LS_plus1Sqrt;
      order_phi_nMinus1_LS=order_phi_nMinus1_LS_plus1Sqrt;
  end
    
  orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

  scalarFluxErrorRMS_plot_handle=figure(iSoln+4);
  
  loglog(gridMeshSize_n,error_phi0_n_FS,'*');
  hold on;
  loglog(gridMeshSize_n,error_phi0_n_LS,'*');
  title({'scalar flux error convergence',[soln ' case']});
  xlabel('mesh size [cm]');
  ylabel('scalar flux error RMS');
  
  orderGuess=round(order_phi_nMinus1_FS(end));
  errorStt=error_phi0_n_FS(end)*refinementRatio^(orderGuess*(nGrids-1));
  firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
  secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
  thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
  fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
  loglog(orderPlotGrid,firstOrder,'r--');
  loglog(orderPlotGrid,secondOrder,'g--');
  loglog(orderPlotGrid,thirdOrder,'b--');
  loglog(orderPlotGrid,fourthOrder,'k--');

  % plot the second set of OoA guide lines
  orderGuess=round(order_phi_nMinus1_LS(end));
  errorStt=error_phi0_n_LS(end)*refinementRatio^(orderGuess*(nGrids-1));
  firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
  secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
  thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
  fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
  loglog(orderPlotGrid,firstOrder,'r--');
  loglog(orderPlotGrid,secondOrder,'g--');
  loglog(orderPlotGrid,thirdOrder,'b--');
  loglog(orderPlotGrid,fourthOrder,'k--');
  
  legend('FS-MoC \phi error','LS-MoC \phi error','1st Order','2nd Order',...
    '3rd Order','4th Order','location','northwest');
  savefig([soln '_withScattering_MC2017']);
  hold off;
  
end
