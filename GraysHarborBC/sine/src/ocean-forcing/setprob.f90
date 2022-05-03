subroutine setprob()

    use setprob_module, only: tide_t, tide_etaprime, tide_m

    implicit none

    real(kind=8) :: tide_eta  ! read in but not used

    integer iunit, i
    character(len=25) fname

    ! if not using tide data...
    return


    iunit = 7
    fname = 'XXX.txt'
    open(unit=iunit, file=fname, status='old')

    read(iunit,*) tide_m
    allocate(tide_t(tide_m), tide_etaprime(tide_m))

    do i=1,tide_m
        read(iunit,*) tide_t(i), tide_eta, tide_etaprime(i)
    enddo

    close(iunit)

    
end subroutine setprob
