subroutine setprob()

    use setprob_module, only: a

    implicit none


    INTEGER :: col,rows,io
    
    rows = 0
    OPEN(UNIT=19, FILE="tidal_signal.txt")
    !count # of rows
    DO
    READ(19,*,iostat=io)
    IF (io/=0) EXIT
    rows = rows + 1
    END DO

    rewind(19)
    !export as a()
    DO col = 1, rows
      READ(19,*) a(col)
    END DO
    close(19)


end subroutine setprob
