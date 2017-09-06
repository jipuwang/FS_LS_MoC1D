%engin

nGrids=8;

assumedSoln='constant';
figureID=11;
[error_phi0_n_FS_constant, order_phi_nMinus1_FS_constant]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_LS_constant, order_phi_nMinus1_LS_constant]=converger_LS(assumedSoln,nGrids,figureID)

assumedSoln='linear';
figureID=12;
[error_phi0_n_FS_linear, order_phi_nMinus1_FS_linear]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_LS_linear, order_phi_nMinus1_LS_linear]=converger_LS(assumedSoln,nGrids,figureID)

assumedSoln='quadratic';
figureID=13;
[error_phi0_n_FS_quadratic, order_phi_nMinus1_FS_quadratic]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_LS_quadratic, order_phi_nMinus1_LS_quadratic]=converger_LS(assumedSoln,nGrids,figureID)

assumedSoln='plus1Sqrt';
figureID=14;
[error_phi0_n_FS_plus1Sqrt, order_phi_nMinus1_FS_plus1Sqrt]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_LS_plus1Sqrt, order_phi_nMinus1_LS_plus1Sqrt]=converger_LS(assumedSoln,nGrids,figureID)

