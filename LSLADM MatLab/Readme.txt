Environment: MATLAB R2022a

The “Example one” folder contains the following files:
1. The file GKDV is used to compute approximate solutions for Example 1 under the integer-order derivative scenario (a=1) and to perform error analysis against the exact solution. It generates plots of the approximate solution, LADM error plots, and LSLADM error plots.
2. The file GKDV2 is used to compute the approximate solution for Example 1 under fractional-order derivative conditions and to perform residual analysis. It generates plots of the approximate solution, LADM residuals, and LSLADM residuals.
3. File GKDV3 is used to calculate the sum of squares of residuals for Example 1 under different fractional-order derivative scenarios. The output results are stored in the function files gkdv1, gkdv08, gkdv05, and gkdv02. This provides support for solving for unknown coefficients using the quasi-Newton method.
4. File quasiNewton is the program for running the quasi-Newton method.

The Example two folder contains the following files:
1. The file ZK is used to compute the approximate solution for Example 1 under the integer-order derivative scenario (a=1) and to perform an error analysis with the exact solution. It generates plots of the approximate solution, LADM error, and LSLADM error.
2. The file ZK2 is used to compute the approximate solution for Example 1 under the fractional-order derivative scenario and to perform a residual analysis. It produces plots of the approximate solution, LADM residuals, and LSLADM residuals.
3. The file ZK3 is used to calculate the sum of squares of residuals for Example 1 under different fractional-order derivative scenarios. The output results are stored in the function files zk1, zk08, zk05, and zk02. This provides support for solving for unknown coefficients using the quasi-Newton method.
4. The file quasiNewton is the program for running the quasi-Newton method.


Usage Instructions (illustrated using Example 1; Example 2 is similar):
1. To compute the solution for integer-order derivatives, run GKDV directly.
2. To compute the solution for fractional-order derivatives, run GKDV2. When calculating the 0.8th-order derivative (a = 0.8), select the ur, resu, and resur values for a = 0.8; comment out the other functions with the same names. The same applies to other fractional-order cases.

Note: GKDV3 and its function files gkdv1, gkdv08, gkdv05, and gkdv02 are backup files. They do not need to be run.
