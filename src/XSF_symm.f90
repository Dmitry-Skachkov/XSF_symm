




     Program XSF_symm                             ! read excitonic wave function from xsf file and make it symmetrical 
      use Module_XSF
      call getarg(1,name)                         ! read name of the file
      Rx0 = 10.5d0                                ! radius of Frenkel exciton
      call read_xsf                               ! read excitonic wave function from file
      call reduce_potential_to_sphere             ! reduce the size of the potential
      call make_symm                              ! make symetric wave function
      call write_xsf                              ! write new wave function
     end program XSF_symm















