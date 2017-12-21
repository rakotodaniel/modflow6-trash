module MpiExchangeGenModule
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize
  
  implicit none
  
  private
  
  ! -- Public functions
  public :: mpi_initialize, mpi_finalize, mpi_append_fname
    
  ! -- Global variables based of the world communicator
  logical, public :: serialrun = .true., writestd = .true.
  character(len=50), public :: partstr
  
  save
  
  contains
  
  subroutine mpi_initialize()
! ******************************************************************************
! MPI initialization to be called at the start of the program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
! ------------------------------------------------------------------------------
    !
    call mpiwrpinit()
    !
    ! -- return
    return
  end subroutine mpi_initialize

  subroutine mpi_finalize()
! ******************************************************************************
! MPI finalization to be called at the end of the program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
! ------------------------------------------------------------------------------
    !
    call mpiwrpfinalize()
    !
    ! -- return
    return
  end subroutine mpi_finalize

  subroutine mpi_append_fname(f)
! ******************************************************************************
! Append the file name with the process rank ID.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    character(len=*), intent(inout) :: f
    ! -- local
! ------------------------------------------------------------------------------
    !
    ! -- Return in case of serial run
    if (serialrun) return
    !
    write(f,'(3a)') trim(f),'.',trim(partstr)
    !
    ! -- return
    return
  end subroutine mpi_append_fname
  
end module MpiExchangeGenModule
