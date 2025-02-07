module talcc_c_bindings

    use, intrinsic :: iso_c_binding, only: c_int

    implicit none

    public

    interface

    subroutine printAffinityBindingC(mpiid, ompid) bind(c, name='talcc_printaffinitymaskstringstdout')
        import :: c_int
        integer(c_int), value :: mpiid
        integer(c_int), value :: ompid
    end subroutine printAffinityBindingC

    subroutine setAffinityBindingC(numcorespernode, numopenmpthreads) bind(c, name='talcc_setprogramaffinityplugin')
        import :: c_int
        integer(c_int), value :: numcorespernode
        integer(c_int), value :: numopenmpthreads
    end subroutine setAffinityBindingC
   
  end interface

end module talcc_c_bindings
