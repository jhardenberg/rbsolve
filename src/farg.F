	integer function mpir_iargc_()
        mpir_iargc_ = iargc()
        return
        end
        subroutine mpir_getarg_( i, s )
        integer       i
        character*(*) s
        call getarg(i,s)
        return
        end


