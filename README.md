# NRG
Numerical renormalization group for IRLM and SAHM
to run the code:
mpif90 -O3 -c mpmodule.f90 mpfuna.f90 mpfunbq.f90 mpfunc.f90 mpfund.f90 mpfune.f90 mpfunf.f90 mpfungq1.f90 second.f90
mpif90 -O3 -o inrg58 inrg58.f90 mpmodule.o mpfuna.o mpfunbq.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfungq1.o second.o -lblas -llapack
./inrg58
