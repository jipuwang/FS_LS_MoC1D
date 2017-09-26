% Need to set the material scattering xs from 0.4 to zero. adn absorption to be
% from 0.5 to 0.9

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
assumedSoln='plus1Sqrt-expMu';
[error_phi0_n_FS_plus1Sqrt, order_phi_nMinus1_FS_plus1Sqrt]=converger_FS(assumedSoln,nGrids);
[error_phi0_n_LS_plus1Sqrt, order_phi_nMinus1_LS_plus1Sqrt]=converger_LS(assumedSoln,nGrids);

% reference solution
%%
% =================
% assumedSoln: constant
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_constant =   1.0e-13 *[...
     0
     0
     0
     0
     0
     0
     0
     0
];

order_phi_nMinus1_FS_constant =0.0*[...
   NaN
   NaN
   NaN
   NaN
   NaN
   NaN
   NaN
];

% =================
% assumedSoln: constant
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_constant =   1.0e-14 *[...
   0.065869074524511
   0.094987349702291
   0.128325645903259
   0.094857496805351
   0.077755252436221
   0.079227462882178
   0.069149829109873
   0.085773261495406
];

order_phi_nMinus1_LS_constant =[...
  -0.528134109823373
  -0.434002227095005
   0.435975819353337
   0.286821663774070
  -0.027060469373325
   0.196274917192361
  -0.310802293614654
];
%%
% =================
% assumedSoln: linear
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_linear =[...
   0.035480928826756
   0.010928673292213
   0.002991448010036
   0.000769892856104
   0.000193990994682
   0.000048595130365
   0.000012154909784
   0.000003039111039
];

order_phi_nMinus1_FS_linear =[...
   1.698925504450935
   1.869202377877907
   1.958114400904940
   1.988668001673360
   1.997106025983069
   1.999272565653576
   1.999817893196999
];

% =================
% assumedSoln: linear
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_linear =   1.0e-13 *[...
   0.051041288575940
   0.153926727010981
   0.194555933918328
   0.197907024795120
   0.224672495448142
   0.231332890713263
   0.294101673071550
   0.411483649226436
];

order_phi_nMinus1_LS_linear =[...
  -1.592507098639821
  -0.337941227658696
  -0.024637839475076
  -0.183000697293464
  -0.042146882039066
  -0.346344589674031
  -0.484520113264912
];

%%
% =================
% assumedSoln: quadratic
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_quadratic =[...
   0.341257092218669
   0.104418458171607
   0.028719233512101
   0.007410170118407
   0.001868591795332
   0.000468180751118
   0.000117110220215
   0.000029281627172
];

order_phi_nMinus1_FS_quadratic =[...
   1.708482268014216
   1.862287610636356
   1.954438676989358
   1.987555225030929
   1.996813912925043
   1.999198635002890
   1.999799354263087
];

% =================
% assumedSoln: quadratic
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_quadratic =[...
   0.002744684878867
   0.000261575875735
   0.000019714240145
   0.000001307909659
   0.000000083070454
   0.000000005213261
   0.000000000326154
   0.000000000020384
];

order_phi_nMinus1_LS_quadratic =[...
   3.391339123730242
   3.729919481854877
   3.913903307022427
   3.976783647943501
   3.994077375866360
   3.998560347708503
   4.000060596171525
];

%%
% =================
% assumedSoln: plus1Sqrt
% refinementRatio: 2
% quad set order: 8

error_phi0_n_FS_plus1Sqrt =[...
   0.005576445250615
   0.001919320485749
   0.000558549438921
   0.000146817122030
   0.000037209391258
   0.000009334933895
   0.000002335786607
   0.000000584075233
];

order_phi_nMinus1_FS_plus1Sqrt =[...
   1.538750126548174
   1.780838742245486
   1.927664756496206
   1.980281533817870
   1.994955077106626
   1.998731327870577
   1.999682363317173
];

% =================
% assumedSoln: plus1Sqrt
% refinementRatio: 2
% quad set order: 8

error_phi0_n_LS_plus1Sqrt =   1.0e-03 *[...
   0.173816204037366
   0.019565401774390
   0.001560753009194
   0.000105189593227
   0.000006708212205
   0.000000421419797
   0.000000026372768
   0.000000001648852
];

order_phi_nMinus1_LS_plus1Sqrt =[...
   3.151185942533062
   3.647990583718791
   3.891178361409248
   3.970919843436341
   3.992598332698723
   3.998137161458585
   3.999514805155529
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
      soln='plus1Sqrt-expMu';
      error_phi0_n_FS=error_phi0_n_FS_plus1Sqrt;
      order_phi_nMinus1_FS=order_phi_nMinus1_FS_plus1Sqrt;
      error_phi0_n_LS=error_phi0_n_LS_plus1Sqrt;
      order_phi_nMinus1_LS=order_phi_nMinus1_LS_plus1Sqrt;
  end
    
  orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

  scalarFluxErrorRMS_plot_handle=figure(iSoln);
  
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
    '3rd Order','4th Order','location','best');
  savefig([soln '_noScattering_MC2017']);
  xlim([0.005 1]);
  hold off;
  
end
