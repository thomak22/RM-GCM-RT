FC = ifort
OPT=-O3 -xHost
#OPT='-Ofast -xHost'
PROF=-g -p
REPORT=-warn noalign
#PARALLEL='-fopenmp'
OTHER_OPTS=-debug extended

FFLAGS = -g $(OPT) $(REPORT) $(PARALLEL) $(OTHER_OPTS) -r8 -132 -traceback

OBJS = 1DRT.o corrk.o ropprmulti_corrk.o radiative_transfer_corrk.o corrk_setup.o roppr1_corrk.o cloud_properties_set_up.o rradiation.o rcalc_radheat.o rradsub.o rsetuprad_simple.o rradtran.o rinterpol.o roppr1.o rtwostr.o radd.o rnewflux1.o rmakeclouds.o ropprmulti.o radiative_transfer_picket_fence.o cbinary.o ciniset.o cset99.o cwrsps.o cinital.o cinigau.o ciniphys.o ciniqs.o cinires.o ciniresij.o cinisi.o inisimprad.o cinistr.o inivarparam.o cgwtlt.o cinibal.o cinisp.o clgndre.o cmatinv.o cqreig.o csetzt.o chessen.o cicamax.o cqrt.o csgetrf.o csgetri.o cilaenv.o csgemm.o csgemv.o csgetf2.o cslaswp.o csswap.o cstrsm.o cstrtri.o cxerbla.o cisamax.o clsame.o csger.o csscal.o cstrmm.o cstrti2.o cstrmv.o

1DRT: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod 1DRT
