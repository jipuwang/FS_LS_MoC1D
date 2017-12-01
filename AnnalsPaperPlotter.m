% Annals paper Plotter
%% Source problem plot
soln='quadratic-expMu';
% refinementRatio to be assigned later
% nGrids: to be assigned later
format long E
%read in data
err_phi0_ij=[
1.00000000000000E+00	1.62529973659790E-02	1.81531494415858E+00	1.62528322242600E-02	1.81535215202695E+00	1.62527259639710E-02	1.81537646100834E+00
5.00000000000000E-01	4.61816255070300E-03	1.94566681176064E+00	4.61799652486100E-03	1.94581502834053E+00	4.61788852190600E-03	1.94591164500739E+00
2.50000000000000E-01	1.19885081622400E-03	1.98493616823559E+00	1.19868456251600E-03	1.98553452228191E+00	1.19857625757000E-03	1.98592181737740E+00
1.25000000000000E-01	3.02858534965000E-04	1.99239813927759E+00	3.02690968946000E-04	1.99488195992667E+00	3.02582379957000E-04	1.99644678853423E+00
6.25000000000000E-02	7.61146428600000E-05	1.98127038253395E+00	7.59416721950000E-05	1.99260246257148E+00	7.58321320290000E-05	1.99910963568855E+00
3.12500000000000E-02	1.92773085770000E-05	1.90583515244349E+00	1.90830173510000E-05	1.96998481834086E+00	1.89697366350000E-05	1.99977727255795E+00
1.56250000000000E-02	5.14437775500000E-06	1.47522905188379E+00	4.87104928600000E-06	1.83454345375622E+00	4.74316636600000E-06	1.99994427783208E+00
7.81250000000000E-03	1.85031068400000E-06	5.40698862966571E-01	1.36574593300000E-06	1.19812423535575E+00	1.18583739200000E-06	1.99998594098755E+00
3.90625000000000E-03	1.27197354900000E-06	9.73483295549306E-02	5.95248874000000E-07	3.31524405670576E-01	2.96462237000000E-07	1.99999607284447E+00
1.95312500000000E-03	1.18897661500000E-06	0.00000000000000E+00	4.73042098000000E-07	0.00000000000000E+00	7.41157610000000E-08	0.00000000000000E+00
  ];

gridMeshSize_n=err_phi0_ij(:,1);
nGrids=size(gridMeshSize_n,1)
refinementRatio=gridMeshSize_n(1)/gridMeshSize_n(2);

err_phi0_no_aer=err_phi0_ij(:,2);
err_phi0_partial_aer=err_phi0_ij(:,4);
err_phi0_complete_aer=err_phi0_ij(:,6);

% the correct OoA
order_phi0_complete_aer=err_phi0_ij(:,7);

% plot for phi0
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

scalarFluxErrorRMS_plot_handle=figure;

loglog(gridMeshSize_n,err_phi0_no_aer,'r*');
hold on;
loglog(gridMeshSize_n,err_phi0_partial_aer,'bx');
loglog(gridMeshSize_n,err_phi0_complete_aer,'ko');

title({'scalar flux error convergence',[soln ' case']});
xlabel('mesh size [cm]');
ylabel('scalar flux error RMS');

orderGuess=round(order_phi0_complete_aer(end-1)); % -1 because I've got an spurious zero.
errorStt=err_phi0_complete_aer(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'r--');
loglog(orderPlotGrid,secondOrder,'g--');
loglog(orderPlotGrid,thirdOrder,'b--');
loglog(orderPlotGrid,fourthOrder,'k--');

legend('no AER','parital AER','complete AER','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');
savefig([soln '_temp_phi0_MoC_AnnalsMMS1D2017']);
hold off;

%% The eigenvalue problem plot
soln='plus1Sqrt-expMu';
% refinementRatio to be assigned later
% nGrids: to be assigned later
format long E
%read in data
err_phi0_ij=[
1.00000000000000E+00	1.62529973659790E-02	1.81531494415858E+00	1.62528322242600E-02	1.81535215202695E+00	1.62527259639710E-02	1.81537646100834E+00
5.00000000000000E-01	4.61816255070300E-03	1.94566681176064E+00	4.61799652486100E-03	1.94581502834053E+00	4.61788852190600E-03	1.94591164500739E+00
2.50000000000000E-01	1.19885081622400E-03	1.98493616823559E+00	1.19868456251600E-03	1.98553452228191E+00	1.19857625757000E-03	1.98592181737740E+00
1.25000000000000E-01	3.02858534965000E-04	1.99239813927759E+00	3.02690968946000E-04	1.99488195992667E+00	3.02582379957000E-04	1.99644678853423E+00
6.25000000000000E-02	7.61146428600000E-05	1.98127038253395E+00	7.59416721950000E-05	1.99260246257148E+00	7.58321320290000E-05	1.99910963568855E+00
3.12500000000000E-02	1.92773085770000E-05	1.90583515244349E+00	1.90830173510000E-05	1.96998481834086E+00	1.89697366350000E-05	1.99977727255795E+00
1.56250000000000E-02	5.14437775500000E-06	1.47522905188379E+00	4.87104928600000E-06	1.83454345375622E+00	4.74316636600000E-06	1.99994427783208E+00
7.81250000000000E-03	1.85031068400000E-06	5.40698862966571E-01	1.36574593300000E-06	1.19812423535575E+00	1.18583739200000E-06	1.99998594098755E+00
3.90625000000000E-03	1.27197354900000E-06	9.73483295549306E-02	5.95248874000000E-07	3.31524405670576E-01	2.96462237000000E-07	1.99999607284447E+00
1.95312500000000E-03	1.18897661500000E-06	0.00000000000000E+00	4.73042098000000E-07	0.00000000000000E+00	7.41157610000000E-08	0.00000000000000E+00
  ];
err_k_ij=[
1.00000000000000E+00	-1.62092553015913E-04	1.37007425740429E+00	-1.61964475436482E-04	1.37188345180368E+00	-1.61885200285816E-04	1.37301034910554E+00
5.00000000000000E-01	-6.27089623947970E-05	1.81909470936221E+00	-6.25808848149220E-05	1.82658005324099E+00	-6.25014146855560E-05	1.83126502788044E+00
2.50000000000000E-01	-1.77716414409000E-05	1.90868281118949E+00	-1.76435638610250E-05	1.93782407450011E+00	-1.75640289423740E-05	1.95644583585593E+00
1.25000000000000E-01	-4.73322079619400E-06	1.81234092882344E+00	-4.60514321654100E-06	1.91683144025509E+00	-4.52559058872200E-06	1.98901105591118E+00
6.25000000000000E-02	-1.34768309911800E-06	1.45025848661800E+00	-1.21960551857600E-06	1.73999199613994E+00	-1.14004835705100E-06	1.99724621961974E+00
3.12500000000000E-02	-4.93192514162000E-07	8.21578543455648E-01	-3.65114934731000E-07	1.27397563462804E+00	-2.85556632784000E-07	1.99931113584052E+00
1.56250000000000E-02	-2.79059422237000E-07	3.07482125236326E-01	-1.50981841474000E-07	6.32135407180467E-01	-7.14232535340000E-08	1.99982766464844E+00
7.81250000000000E-03	-2.25494185546000E-07	8.83396331818082E-02	-9.74166070030000E-08	2.13378304234262E-01	-1.78579464550000E-08	1.99995632083339E+00
3.90625000000000E-03	-2.12100879082000E-07	2.29576348649451E-02	-8.40232989850000E-08	5.86704979797296E-02	-4.46462178300000E-09	1.99999318403751E+00
1.95312500000000E-03	-2.08752426234000E-07	0.00000000000000E+00	-8.06748468030000E-08	0.00000000000000E+00	-1.11616071900000E-09	0.00000000000000E+00
  ];

gridMeshSize_n=err_phi0_ij(:,1);
nGrids=size(gridMeshSize_n,1)
refinementRatio=gridMeshSize_n(1)/gridMeshSize_n(2);

err_phi0_no_aer=err_phi0_ij(:,2);
err_k_no_aer=err_k_ij(:,2);
err_phi0_partial_aer=err_phi0_ij(:,4);
err_k_partial_aer=err_k_ij(:,4);
err_phi0_complete_aer=err_phi0_ij(:,6);
err_k_complete_aer=err_k_ij(:,6);

% the correct OoA
order_phi0_complete_aer=err_phi0_ij(:,7);
order_k_complete_aer=err_k_ij(:,7);

% plot for phi0
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

scalarFluxErrorRMS_plot_handle=figure;

loglog(gridMeshSize_n,err_phi0_no_aer,'r*');
hold on;
loglog(gridMeshSize_n,err_phi0_partial_aer,'bx');
loglog(gridMeshSize_n,err_phi0_complete_aer,'ko');

title({'scalar flux error convergence',[soln ' case']});
xlabel('mesh size [cm]');
ylabel('scalar flux error RMS');

orderGuess=round(order_phi0_complete_aer(end-1)); % -1 because I've got an spurious zero.
errorStt=err_phi0_complete_aer(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'r--');
loglog(orderPlotGrid,secondOrder,'g--');
loglog(orderPlotGrid,thirdOrder,'b--');
loglog(orderPlotGrid,fourthOrder,'k--');

legend('no AER','parital AER','complete AER','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');

% set(get(gca,'xlabel'),'FontName','Times New Roman');
% set(get(gca,'ylabel'),'FontName','Times New Roman');
% set(get(gca,'title'),'FontName','Times New Roman');
% set(findobj(gcf, 'Type', 'Legend'),'FontName','Times New Roman');

savefig([soln '_temp_phi0_MoCEign_AnnalsMMS1D2017']);
hold off;

% plot for k
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

k_Error_plot_handle=figure;

loglog(gridMeshSize_n,err_k_no_aer,'r*');
hold on;
loglog(gridMeshSize_n,err_k_partial_aer,'bx');
loglog(gridMeshSize_n,err_k_complete_aer,'ko');

title({'k error convergence',[soln ' case']});
xlabel('mesh size [cm]');
ylabel('k error');

orderGuess=round(order_k_complete_aer(end-1)); % -1 because I've got an spurious zero.
errorStt=err_k_complete_aer(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'r--');
loglog(orderPlotGrid,secondOrder,'g--');
loglog(orderPlotGrid,thirdOrder,'b--');
loglog(orderPlotGrid,fourthOrder,'k--');

legend('no AER','parital AER','complete AER','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');

% set(get(gca,'xlabel'),'FontName','Times New Roman');
% set(get(gca,'ylabel'),'FontName','Times New Roman');
% set(get(gca,'title'),'FontName','Times New Roman');
% set(findobj(gcf, 'Type', 'Legend'),'FontName','Times New Roman');

savefig([soln '_temp_k_MoCEign_AnnalsMMS1D2017']);
hold off;

%%