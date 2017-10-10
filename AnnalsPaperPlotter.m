% Annals paper Plotter
%% Problem discription
soln='quadratic-expMu';
refinementRatio=4;
% nGrids: to be assigned later
%% The eigenvalue plot
format long E
%read in data
err_phi0_ij=aa;
err_k_ij=bb;

gridMeshSize_n=err_phi0_ij(:,1);
nGrids=size(gridMeshSize_n,1)

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

loglog(gridMeshSize_n,err_phi0_no_aer,'*');
hold on;
loglog(gridMeshSize_n,err_phi0_partial_aer,'*');
loglog(gridMeshSize_n,err_phi0_complete_aer,'*');

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
savefig([soln '_phi0_MoCEign_AnnalsMMS1D2017']);
hold off;

% plot for k
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

k_Error_plot_handle=figure;

loglog(gridMeshSize_n,err_k_no_aer,'*');
hold on;
loglog(gridMeshSize_n,err_k_partial_aer,'*');
loglog(gridMeshSize_n,err_k_complete_aer,'*');

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
savefig([soln '_k_MoCEign_AnnalsMMS1D2017']);
hold off;


%%