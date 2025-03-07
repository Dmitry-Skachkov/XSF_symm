default: xsf_symm

FC = ifort
FFLAGS = 

obj/Qshep.o: Qshep.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -c Qshep.f90           -o obj/Qshep.o        -module obj
obj/Module_XSF.o: Module_XSF.f90
	$(FC) $(FFLAGS) -c Module_XSF.f90      -o obj/Module_XSF.o   -module obj
obj/XSF_symm.o: XSF_symm.f90
	$(FC) $(FFLAGS) -c XSF_symm.f90        -o obj/XSF_symm.o     -module obj


OBJS =    obj/Qshep.o          \
          obj/Module_XSF.o     \
          obj/XSF_symm.o


xsf_symm: $(OBJS)
	$(FC) $(LNK) -o ~/bin/xsf_symm $(OBJS) -module obj


clean:
	rm -rf ~/bin/xsf_symm obj/*.o obj/*.mod

