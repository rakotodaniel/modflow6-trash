module MpiWrapper
  
  use ConstantsModule, only: LENORIGIN, LENVARNAME
  use KindModule, only: DP, I4B    
  
  implicit none
 
#ifdef MPI_PARALLEL
  include 'mpif.h'
#endif
  
  private
  
  integer, parameter :: mpiwrplblk_size = 64
  integer, parameter :: mpiwrprbuf_size = 256
#ifdef MPI_PARALLEL
  character(len=MPI_MAX_PROCESSOR_NAME) :: mpiwrpmyname
  integer, dimension(MPI_STATUS_SIZE,mpiwrplblk_size) :: mpiwrpstats
#endif  
  integer :: mpiwrpcomm_world, mpiwrpcomm_null
  integer :: mpiwrpnrproc, mpiwrpmyrank
  integer, dimension(mpiwrprbuf_size) :: mpiwrprbufi
  double precision, dimension(mpiwrprbuf_size) :: mpiwrprbufd
  integer, dimension(1000) :: rreq, sreq
  integer :: lenbuf
  
  type ColMemoryType
    character(len=LENVARNAME)  :: name        !name of the array
    character(len=LENORIGIN)   :: origin      !name of origin
    integer(I4B)               :: memitype    !integer type
    logical                    :: logicalsclr !logical
    integer(I4B)               :: intsclr     !integer
    real(DP)                   :: dblsclr     !double
    integer(I4B), dimension(3) :: aint1d      !1d integer array
  end type ColMemoryType
  
  public :: mpiwrpcomm_size, mpiwrpcomm_rank, mpiwrpinit, mpiwrpbarrier
  public :: mpiwrpfinalize, mpiwrpcommworld
  public :: mpiwrpnrproc, mpiwrpmyrank
  public :: mpiwrpisend, mpiwrpallsum
  public :: mpiwrpallgather, mpiwrpallgatherv, mpiwrpcommgroup
  public :: mpiwrpgroupincl, mpiwrpcommcreate
  public :: mpiwrpgrouprank
  public :: mpiwrpcolstruct
  public :: ColMemoryType
  public :: mpiwrptypefree
  
  save
  
  interface mpiwrpisend
    module procedure mpiwrpisendd
  end interface
  
  interface mpiwrpallsum
    module procedure mpiwrpallsumr, mpiwrpallsumd
  end interface
  
  interface mpiwrpallgather
    module procedure mpiwrpallgatheri
  end interface 
  
  interface mpiwrpallgatherv
    module procedure mpiwrpallgathervi, mpiwrpallgathervd,                      &
                     mpiwrpallgathervcol 
  end interface
  
  contains
  
  integer function mpiwrpcomm_size(comm)
! ******************************************************************************
! Wrapper subroutine around MPI_COMM_SIZE te determine number of
! processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm ! (I) communicator
    ! -- local
    integer :: ierror, size
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_comm_size(comm, size, ierror)
    mpiwrpcomm_size = size
#else
    mpiwrpcomm_size = 1
#endif
    ! -- return
    return
  end function mpiwrpcomm_size
  
  integer function mpiwrpcomm_rank(comm)
! ******************************************************************************
! Wrapper function around MPI_COMM_RANK to determine own
! process rank ID.
! @return rank ID of my own process.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm ! (I) communicator
    ! -- local
    integer :: ierror, rank
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
   call mpi_comm_rank( comm, rank, ierror )
   mpiwrpcomm_rank = rank
#else
   mpiwrpcomm_rank = 0
#endif
    ! -- return
    return
  end function mpiwrpcomm_rank

   subroutine mpiwrpinit()
! ******************************************************************************
! Wrapper subroutine around MPI_INIT for initialization.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierror, i, required, provided
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL    
    required = MPI_THREAD_FUNNELED
    call mpi_init_thread(required, provided, ierror)
    
    if (required /= provided) then
      if (mpiwrpmyrank == 0) then
        write(*,*) 'Warning, could not guarantee thread safety.'
      end if
    end if   
#endif
    ! -- return
    return
  end subroutine mpiwrpinit

  subroutine mpiwrpfinalize()
! ******************************************************************************
! Wrapper subroutine around MPI_FINALIZE termination.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_finalize(ierror)
#endif
    ! -- return
    return
  end subroutine mpiwrpfinalize
  
  subroutine mpiwrpbarrier(comm)
! ******************************************************************************
! Wrapper subroutine around MPI_BARRIER barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_barrier(comm, ierror)
#endif
    ! -- return
    return
  end subroutine mpiwrpbarrier
  
  integer function mpiwrpcommworld()
! ******************************************************************************
! Wrapper subroutine around MPI_BARRIER barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
  
#ifdef MPI_PARALLEL
    mpiwrpcommworld = MPI_COMM_WORLD
#else
    mpiwrpcommworld = 0
#endif
    ! -- return
    return
  end function mpiwrpcommworld

  subroutine mpiwrpisendd(dbuf, count, dest, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking send of a double array to the
! process with rankID dest.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    double precision, dimension(*), intent(in) :: dbuf  ! double buffer to be sent
    integer, intent(in)                        :: count ! number of integers to be sent
    integer, intent(in)                        :: dest  ! rank id of destination process
    integer, intent(in)                        :: tag   ! message tag
    integer, intent(in)                        :: comm  ! communicator
    integer, intent(out)                       :: req   ! request handle
    ! -- local
    integer :: ierror      
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_isend(dbuf, count, mpi_double_precision, dest, tag,               &
                   comm, req, ierror)
#else
    req = -1
    call mpiwrperror(comm, 'mpiwrpisendd', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpisendd
  
  subroutine mpiwrpirecvd(dbuf, count, source, tag, comm, req)
! ******************************************************************************
! Begins a non-blocking receive of an double array from
! the process with rank source.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    double precision, dimension(*), intent(out) :: dbuf   ! (I) double buffer
    integer, intent(in)                         :: count  ! (I) number of integers to be received
    integer, intent(in)                         :: source ! (I) rank of source process
    integer, intent(in)                         :: tag    ! (I) message tag
    integer, intent(in)                         :: comm   ! (I) communicator
    integer, intent(out)                        :: req    ! (O) request handle
    ! -- local
    integer  :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
      call mpi_irecv(dbuf, count, mpi_double_precision, source, tag,           &
                     comm, req, ierror)
#else
      call mpiwrperror(comm, 'mpiwrpirecvd', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpirecvd

  subroutine mpiwrpallsumr(comm, gsbuf, grbuf, n)
! ******************************************************************************
! Communicates and adds double value(s) across all processes.
! To add a scalar value, simply use gsbuf with size 1.
! For arrays, each index is minimized separately (across processes).
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm  ! Communicator
    integer, intent(in) :: n     ! Number of elements in data array.
    real, dimension(n)  :: gsbuf ! Array with value(s) to be minimized.
    real, dimension(n)  :: grbuf ! Temporary receive buffer.
    ! -- local
    integer  :: ierror, i
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allreduce(gsbuf,                                                  &
                       grbuf,                                                  &
                       n,                                                      &
                       mpi_real,                                               &
                       mpi_sum,                                                &
                       comm,                                                   &
                       ierror)
    !
    do i = 1, n
      gsbuf(i) = grbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallsumd', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallsumr
  
  subroutine mpiwrpallsumd(comm, gsbuf, grbuf, n)
! ******************************************************************************
! Communicates and adds double value(s) across all processes.
! To add a scalar value, simply use gsbuf with size 1.
! For arrays, each index is minimized separately (across processes).
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)            :: comm  ! Communicator
    integer, intent(in)            :: n     ! Number of elements in data array.
    double precision, dimension(n) :: gsbuf ! Array with value(s) to be minimized.
    double precision, dimension(n) :: grbuf ! Temporary receive buffer.
    ! -- local
    integer  :: ierror, i
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allreduce(gsbuf,                                                  &
                       grbuf,                                                  &
                       n,                                                      &
                       mpi_double_precision,                                   &
                       mpi_sum,                                                &
                       comm,                                                   &
                       ierror)
    !
    do i = 1, n
      gsbuf(i) = grbuf(i)
    end do
    !
#else
    call mpiwrperror(comm, 'mpiwrpallsumd', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpallsumd
  
  subroutine mpiwrperror(comm, subname, message)
! ******************************************************************************
! Prints an error message and aborts the current program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in)          :: comm    ! (I) communicator
    character(len=*), intent(in) :: subname ! (I) the name of the subroutine in which the error occured
    character(len=*), intent(in) :: message ! (I) the error message to be printed
    ! -- local
    integer :: ierror, i
    integer, dimension(1) :: idummy(1)
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    write (6,*)
    write (6,*) '*** fatal runtime error in subroutine '//subname
    write (6,*)
    write (6,*) '    '//message
    write (6,*)
    write (6,*) '*** trying to terminate all processes...'
    write (6,*)
    !
    call mpi_abort(comm, 1, ierror)
#else
    write (6,*)
    write (6,*) '*** fatal runtime error in subroutine '//subname
    write (6,*)
    write (6,*) '    '//message
    write (6,*)
    write (6,*) '*** aborting program...'
    write (6,*)
    !
    ! --   try to force a core dump
    i         = -666666666
    idummy(i) =  666666666
    !
    stop
    !
#endif
    ! -- return
    return
  end subroutine mpiwrperror

  subroutine mpiwrpallgatheri(comm, gsbuf, gscnt, grbuf, grcnt)
! ******************************************************************************
! Gathers integer values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer :: grcnt
    integer, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierror, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allgather(gsbuf, gscnt, mpi_integer,                               &
                       grbuf, grcnt, mpi_integer,                               &
                       comm, ierror)
#endif
    ! -- return
    return
  end subroutine mpiwrpallgatheri

  subroutine mpiwrpallgathervcol(comm, gsbuf, gscnt, gstype, grbuf, grcnt,      &
                                 grtype, offsets)
! ******************************************************************************
! Gathers integer values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt, gstype, grtype
    integer, dimension(*) :: grcnt
    type(ColMemoryType), dimension(*) :: gsbuf
    type(ColMemoryType), dimension(*) :: grbuf
    integer, dimension(*) :: offsets

    ! -- local
    integer :: ierror, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allgatherv(gsbuf, gscnt, gstype,                                  &
                        grbuf, grcnt, offsets, grtype,                         &
                        comm, ierror)
    !
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervcol
  
  subroutine mpiwrpallgathervi(comm, gsbuf, gscnt, grbuf, grcnt, offsets)
! ******************************************************************************
! Gathers integer values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer, dimension(*) :: grcnt, offsets
    integer, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierror, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allgatherv(gsbuf, gscnt, mpi_integer,                             &
                        grbuf, grcnt, offsets, mpi_integer,                    &
                        comm, ierror)
#else
    if ( grcnt(1) < gscnt ) then
      call mpiwrperror(comm, 'mpiwrpallgathervi',                              &
                        'receive buffer too small')
    endif
    !
    k = offsets(1)
    !
    do i = 1, gscnt
      grbuf(k+i) = gsbuf(i)
    enddo
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervi
  
  subroutine mpiwrpallgathervd(comm, gsbuf, gscnt, grbuf, grcnt, offsets)
! ******************************************************************************
! Gathers double values into specified memory
! locations from a group of processes and broadcasts the
! gathered data to all processes.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer :: comm, gscnt
    integer, dimension(*) :: grcnt, offsets
    double precision, dimension(*) :: gsbuf, grbuf
    ! -- local
    integer :: ierror, i, k
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_allgatherv(gsbuf, gscnt, mpi_double_precision,                    &
                        grbuf, grcnt, offsets, mpi_double_precision,           &
                        comm, ierror)
#else
    if ( grcnt(1) < gscnt ) then
      call mpiwrperror(comm, 'mpiwrpallgathervd',                              &
                       'receive buffer too small')
    endif
    !
    k = offsets(1)
    !
    do i = 1, gscnt
      grbuf(k+i) = gsbuf(i)
    enddo
#endif
    ! -- return
    return
  end subroutine mpiwrpallgathervd
  
  subroutine mpiwrpcommgroup(comm, group)
! ******************************************************************************
! Accesses the group associated with given communicator 
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: comm
    integer, intent(out) :: group
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_comm_group(comm, group, ierror)
#else
    call mpiwrperror(comm, 'mpiwrpcommgroup', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcommgroup

  subroutine mpiwrpgroupincl(old_group, n, rnks, new_group)
! ******************************************************************************
! Produces a group by reordering an existing group and taking only listed 
! members.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: old_group
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: rnks
    integer, intent(out) :: new_group
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_group_incl(old_group, n, rnks, new_group, ierror)
#else
    call mpiwrperror(comm, 'mpiwrpgroupincl', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpgroupincl
  
  subroutine mpiwrpcommcreate(old_comm, group, new_comm)
! ******************************************************************************
! Creates a new communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: old_comm
    integer, intent(in) :: group
    integer, intent(out) :: new_comm
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_comm_create(old_comm, group, new_comm, ierror)
#else
    call mpiwrperror(comm, 'mpiwrpcommcreate', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcommcreate
  
  subroutine mpiwrpgrouprank(group, rank)
! ******************************************************************************
! Returns the rank of this process in the given group.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: group
    integer, intent(out) :: rank
    ! -- local
    integer :: ierror
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
    call mpi_group_rank(group, rank, ierror)
#else
    call mpiwrperror(comm, 'mpiwrpgrouprank', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpgrouprank

  subroutine mpiwrpcolstruct(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: LENORIGIN, LENVARNAME
    use KindModule, only: DP, I4B    
    ! -- dummy
    integer, intent(out) :: newtype
    ! -- local
    type(ColMemoryType) :: cmtdum
    
    integer, dimension(7) :: types, blocklengths
    integer(KIND=MPI_ADDRESS_KIND), dimension(7) :: displacements
    integer(KIND=MPI_ADDRESS_KIND) :: base
    
    integer :: i, ierr 
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL

    ! -- set-up derived MPI data type 
    types(1) = MPI_CHARACTER
    types(2) = MPI_CHARACTER
    types(3) = MPI_INTEGER
    types(4) = MPI_INTEGER
    types(5) = MPI_REAL
    types(6) = MPI_DOUBLE
    types(7) = MPI_INTEGER
    blocklengths(1)  = LENVARNAME
    blocklengths(2)  = LENORIGIN
    blocklengths(3)  = 1
    blocklengths(4)  = 1
    blocklengths(5)  = 1
    blocklengths(6)  = 1
    blocklengths(7)  = 3
    call MPI_GET_ADDRESS(cmtdum%name,       displacements(1), ierr)
    call MPI_GET_ADDRESS(cmtdum%origin,     displacements(2), ierr)
    call MPI_GET_ADDRESS(cmtdum%memitype,   displacements(3), ierr)
    call MPI_GET_ADDRESS(cmtdum%logicalsclr,displacements(4), ierr)
    call MPI_GET_ADDRESS(cmtdum%intsclr,    displacements(5), ierr)
    call MPI_GET_ADDRESS(cmtdum%dblsclr,    displacements(6), ierr)
    call MPI_GET_ADDRESS(cmtdum%aint1d,     displacements(7), ierr)
    base = displacements(1)
    do i = 1, 7
      displacements(i) = displacements(i) - base
    enddo
    ! -- create and commit datatype
    call MPI_TYPE_CREATE_STRUCT(7, blocklengths, displacements, types, newtype, ierr)
    call MPI_TYPE_COMMIT(newtype, ierr)
#else
    call mpiwrperror(comm, 'mpiwrpstruct', 'invalid operation')
#endif
    ! -- return
    return
  end subroutine mpiwrpcolstruct

  subroutine mpiwrptypefree(newtype)
! ******************************************************************************
! Create a new MPI data type.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    integer, intent(in) :: newtype
    ! -- local
    integer :: ierr
! ------------------------------------------------------------------------------
#ifdef MPI_PARALLEL
  
   call MPI_Type_free(newtype, ierr)
#endif
    ! -- return
    return
  end subroutine mpiwrptypefree
  
end module MpiWrapper