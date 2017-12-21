module MpiExchangeModule
  use MpiExchangeGenModule, only: serialrun, writestd, partstr
  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, LENORIGIN,          &
                             LINELENGTH
  use gwfModule, only: gwfModelType
   use ListModule, only: ListType
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize, mpiwrpnrproc, mpiwrpmyrank,&
                        mpiwrpbarrier, mpiwrpcommworld, mpiwrpcomm_size,       &
                        mpiwrpcomm_rank, mpiwrpallgatherv
  
  implicit none
  
  private
  
  ! -- Public types
  public MpiExchangeType
  ! -- Public variables
  public MpiWorld
  ! -- Public functions
  public mpi_initialize_world
    
  ! -- Local types
  type :: MpiGwfCommInt
    integer(I4B), pointer                                :: xprnk        => null() ! neighboring rank
    character(len=LENPACKAGENAME), dimension(:), pointer :: name         => null() ! variable name
    character(len=LENORIGIN), dimension(:), pointer      :: origin       => null() ! variable origin
    integer(I4B), dimension(:), pointer                  :: xpbuf        => null() ! buffer pointer
    type(gwfModelType), pointer                          :: sndgwfmodel1 => null() ! send GWF model
    type(gwfModelType), pointer                          :: recgwfmodel2 => null() ! receive GWF model
  end type MpiGwfCommInt
      
  type :: MpiExchangeType
    character(len=LENPACKAGENAME)                          :: name          ! name (origin)
    integer                                                :: gnmodel = 0   ! number of global models
    character(len=LENMODELNAME), dimension(:), allocatable :: gmodelnames   ! global model names
    integer                                                :: gnidsoln = 0  ! number of global solution ids
    integer(I4B), dimension(:), allocatable                :: gidsoln       ! global solution ids
    integer                                                :: gnsub = 0     ! number of global subdomains
    integer(I4B), dimension(:), allocatable                :: gsubs         ! global subdomains
    integer                                                :: lnmodel = 0   ! number of llobal models
    character(len=LENMODELNAME), dimension(:), allocatable :: lmodelnames   ! local model names
    integer                                                :: lnsub = 0     ! number of local subdomains
    integer(I4B), dimension(:), allocatable                :: lsubs         ! model local subdomains
    type(ListType), pointer                                :: lmodellist    ! local model list
    type(ListType), pointer                                :: lexchangelist ! local exchanges list 
    integer(I4B), pointer                      :: comm      => null() ! MPI communicator
    integer(I4B), pointer                      :: nrproc    => null() ! number of MPI process for this communicator
    integer(I4B), pointer                      :: myrank    => null() ! MPI rank in for this communicator
    integer(I4B), pointer                      :: myproc    => null() ! MPI proc in for this communicator
    integer(I4B), dimension(:), pointer        :: sbufi     => null() ! MPI send buffer (integer)
    real(DP), dimension(:), pointer            :: sbufd     => null() ! MPI send buffer (double)
    integer(I4B), dimension(:), pointer        :: topolia   => null() ! Topology IA array of size nrproc+1
    integer(I4B), dimension(:), pointer        :: topolja   => null() ! Topology JA array
    integer(I4B), pointer                      :: nrxp      => null() ! Number of exchange partners
    type(MpiGwfCommInt), dimension(:), pointer :: xp     => null() ! MPI data structure point-to-point communication
    character(len=50), pointer                 :: nrprocstr => null() ! Number of processes string
    integer(I4B), pointer                      :: npdigits  => null() ! Number of digits for nrproc
    character(len=50), pointer                 :: partstr   => null() ! Partition string
  contains
    procedure :: mpi_barrier
    procedure :: mpi_create_output_str
    procedure :: mpi_append_fname
    procedure :: mpi_is_iproc
    procedure :: mpi_addmodel
    procedure :: mpi_addsub
    procedure :: mpi_addidsoln
    procedure :: mpi_local_exchange_init
    procedure :: mpi_debug
    procedure :: mpi_da
  end type MpiExchangeType
  
  ! -- World communicator
  type(MpiExchangeType), pointer :: MpiWorld => null()
    
  save
  
  contains
  
  subroutine mpi_initialize_world()
! ******************************************************************************
! MPI initialization for the world communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    ! -- local
    character(len=LENORIGIN) :: origin
! ------------------------------------------------------------------------------
    !
    ! -- Allocate MpiWorld object
    allocate(MpiWorld)
    !
    ! -- Set name
    MpiWorld%name = 'MPI_WORLD'
    !
    ! -- Allocate scalars
    origin = MpiWorld%name
    call mem_allocate(MpiWorld%comm, 'COMM', origin)
    call mem_allocate(MpiWorld%nrproc, 'NRPROC', origin)
    call mem_allocate(MpiWorld%myrank, 'MYRANK', origin)
    call mem_allocate(MpiWorld%myproc, 'MYPROC', origin)
    !
    MpiWorld%comm   = mpiwrpcommworld()
    MpiWorld%nrproc = mpiwrpcomm_size(MpiWorld%comm)
    MpiWorld%myrank = mpiwrpcomm_rank(MpiWorld%comm)
    MpiWorld%myproc = MpiWorld%myrank + 1
    !
    if (MpiWorld%nrproc == 1) then
      serialrun = .true.
    else
      serialrun = .false.
    endif
    ! @@@@@@ DEBUG
    serialrun = .false.
    if (MpiWorld%myrank == 0) then
      writestd = .true.
    else
      writestd = .false.
    endif
    !
    ! -- output strings
    call MpiWorld%mpi_create_output_str()
    partstr = MpiWorld%partstr
    !
    ! -- return
    return
  end subroutine mpi_initialize_world
  
 subroutine mpi_create_comm(this)
! ******************************************************************************
! MPI communicator construction.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    integer, dimension(:), allocatable :: rnks
! ------------------------------------------------------------------------------
    !
    allocate(this%comm)
    ! 
    this%comm = mpiwrpcommworld() ! TODO: support different groups
    !
    ! -- return
    return
 end subroutine mpi_create_comm
   
  subroutine mpi_barrier(this)
! ******************************************************************************
! MPI barrier.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    !
    call mpiwrpbarrier(this%comm)
    !
    ! -- return
    return
  end subroutine mpi_barrier
  
  subroutine mpi_create_output_str(this)
! ******************************************************************************
! Create several strings for output.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    character(len=LENORIGIN) :: origin
    character(len=20) :: fmt  
! ------------------------------------------------------------------------------
    ! -- Allocate
    origin = this%name
    allocate(this%nrprocstr)
    call mem_allocate(this%npdigits, 'NPDIGITS', origin)
    allocate(this%partstr)
    !
    write(this%nrprocstr,*) this%nrproc
    this%nrprocstr = adjustl(this%nrprocstr)
    this%npdigits = len_trim(this%nrprocstr)
    write(fmt,*) this%npdigits
    fmt = adjustl(fmt)
    fmt = '(a,i'//trim(fmt)//'.'//trim(fmt)//')'
    write(this%partstr,fmt) 'p',this%myrank
    !    
    ! -- return
    return
  end subroutine mpi_create_output_str
  
  subroutine mpi_append_fname(this,f)
! ******************************************************************************
! Append the file name with the process rank ID.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*), intent(inout) :: f
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) return
    !
    write(f,'(3a)') trim(f),'.',trim(this%partstr)
    !
    ! -- return
    return
  end subroutine mpi_append_fname
  
  subroutine mpi_debug(this)
! ******************************************************************************
! Debug subroutine.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
! ------------------------------------------------------------------------------
    if (this%myrank == 0) then
      write(*,*) 'Press a key...' 
      pause
    end if
    call mpi_barrier(this)
    ! -- return
    return
  end subroutine mpi_debug
    
  subroutine mpi_local_exchange_init(this)
! ******************************************************************************
! This subroutine initializes the local exchange data structure for the
  ! MPI_COMM_WORLD communicator.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use NumericalModelModule, only: NumericalModelType,                        &
                                    GetNumericalModelFromList
    use NumericalExchangeModule, only: NumericalExchangeType,                  &
                                    GetNumericalExchangeFromList
    use BaseModelModule, only: GetBaseModelFromList
    use ArrayHandlersModule, only: ifind
    use GwfModule, only: GwfModelType
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    class(NumericalModelType), pointer :: mp
    class(NumericalExchangeType), pointer :: cp
    integer(I4B) :: im, ic, i, isub1, isub2
    type(GwfModelType), pointer :: m1, m2
! ------------------------------------------------------------------------------
    !
    ! -- loop over the models within this solution
    do im=1,this%lmodellist%Count()
      mp => GetNumericalModelFromList(this%lmodellist, im)
      ! -- add model to local MPI model list
      call this%mpi_addmodel(2, mp%name)
      ! -- find the subdomain for this modell
      i = ifind(MpiWorld%gmodelnames, mp%name)
      isub1 = MpiWorld%gsubs(i)
      ! -- add subdomain to local MPI subdomain list
      call this%mpi_addsub(2, isub1)
    enddo
    !
    ! -- loop over exchanges
    do ic=1,this%lexchangelist%Count()
      cp => GetNumericalExchangeFromList(this%lexchangelist, ic)
      call cp%get_m1m2(m1, m2)
      ! -- Add to interface
      if (cp%m2_ishalo) then
        i = ifind(MpiWorld%gmodelnames, m1%name)
        isub1 = MpiWorld%gsubs(i)
        i = ifind(MpiWorld%gmodelnames, m2%name)
        isub2 = MpiWorld%gsubs(i)
      end if
      write(*,*) '@@@ m1:',m1%name
      write(*,*) '@@@ m2:',m2%name, cp%m2_ishalo
    enddo    
    !
    ! -- return
    return
  end subroutine mpi_local_exchange_init
  
  subroutine mpi_is_iproc(this, isub, add)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    ! -- dummy
    class(MpiExchangeType) :: this
    integer, intent(in)  :: isub
    logical, intent(out) :: add
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) then
      add = .true.
      return
    endif
    !
    if (isub == this%myproc) then
      add = .true.
    else
      add = .false.
    end if  
    !
    ! -- return
    return
  end subroutine mpi_is_iproc

  subroutine mpi_addmodel(this, iopt, mname)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: iopt
    character(len=*), intent(in) :: mname
    ! -- local
    integer(I4B) :: n
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    endif
    !
    select case(iopt)
      ! store global
      case(1)
        call ExpandArray(this%gmodelnames)
        n = this%gnmodel
        n = n + 1
        this%gmodelnames(n) = mname
        this%gnmodel = n
      ! store local
      case(2)
        call ExpandArray(this%lmodelnames)
        n = this%lnmodel
        n = n + 1
        this%lmodelnames(n) = mname
        this%lnmodel = n
      case default
        write(errmsg,'(a)') 'Program error in mpi_addmodel.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addmodel
  
  subroutine mpi_addsub(this, iopt, isub)
! ******************************************************************************
! Add a model to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: iopt
    integer(I4B), intent(in) :: isub
    ! -- local
    integer(I4B) :: n
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    endif
    !
    select case(iopt)
      ! store global
      case(1)
        call ExpandArray(this%gsubs)
        n = this%gnsub
        n = n + 1
        this%gsubs(n) = isub
        this%gnsub = n
      ! store local
      case(2)
        call ExpandArray(this%lsubs)
        n = this%lnsub
        n = n + 1
        this%lsubs(n) = isub
        this%lnsub = n
      case default
        write(errmsg,'(a)') 'Program error in mpi_addsub.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addsub
  
  subroutine mpi_addidsoln(this, iopt, idsoln)
! ******************************************************************************
! Add a solution ID to my subdomain.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ArrayHandlersModule, only: ExpandArray
    use SimModule, only: store_error, ustop
    ! -- dummy
    class(MpiExchangeType) :: this
    integer(I4B), intent(in) :: iopt
    integer(I4B), intent(in) :: idsoln
    ! -- local
    integer(I4B) :: n
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    if (serialrun) then
      return
    endif
    !
    select case(iopt)
      ! store global
      case(1)
        call ExpandArray(this%gidsoln)
        n = this%gnidsoln
        n = n + 1
        this%gidsoln(n) = idsoln
        this%gnidsoln = n
      case default
        write(errmsg,'(a)') 'Program error in mpi_addidsoln.'
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_addidsoln
  
  subroutine mpi_da(this)
! ******************************************************************************
! mpi_da -- deallocate
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(MpiExchangeType) :: this
! ------------------------------------------------------------------------------
    !
    call mem_deallocate(this%comm)
    call mem_deallocate(this%nrproc)
    call mem_deallocate(this%myrank)
    call mem_deallocate(this%myproc)
    call mem_deallocate(this%sbufi)
    call mem_deallocate(this%sbufd)
    call mem_deallocate(this%topolia)
    call mem_deallocate(this%topolja)
    call mem_deallocate(this%nrxp)
    deallocate(this%nrprocstr)
    call mem_deallocate(this%npdigits)
    deallocate(this%partstr)
    !
    ! -- return
    return
  end subroutine mpi_da

  !  subroutine mpi_get_dis_world(nlm, lmidx, lmdistype, ngm, gmdistype)
!! ******************************************************************************
!! This subroutine gathers DIS information for all models.
!! ******************************************************************************
!!
!!    SPECIFICATIONS:
!! ------------------------------------------------------------------------------
!    ! -- modules
!    ! -- dummy
!    integer, intent(in)  :: nlm
!    integer, dimension(nlm), intent(in) :: lmidx
!    integer, dimension(nlm), intent(in) :: lmdistype
!    integer, intent(in)  :: ngm
!    integer, dimension(ngm), intent(out) :: gmdistype
!    ! -- local
!    integer :: iproc, jproc, scnt, i, j
!    integer, dimension(:), allocatable :: sbuf, rbuf, rcnt, offsets
!! ------------------------------------------------------------------------------
!    !
!    ! -- allocate work arrays
!    allocate(sbuf(ngm), rbuf(world_nrproc*ngm))
!    allocate(rcnt(world_nrproc), offsets(world_nrproc))
!    !
!    scnt = ngm
!    sbuf = -1
!    do i = 1, nlm
!      j = lmidx(i)
!      sbuf(j) = lmdistype(i)
!    enddo
!    !
!    ! -- offsets
!    offsets(1) = 0
!    do iproc = 1, world_nrproc-1
!      offsets(iproc+1) = offsets(iproc) + rcnt(iproc)
!    end do
!    !
!    ! -- MPI all-to-all
!    call mpiwrpallgatherv(world_comm, sbuf, scnt, rbuf, rcnt, offsets)
!    !
!    gmdistype = -1
!    do iproc = 1, world_nrproc
!      i = offsets(iproc)
!      do jproc = 1, world_nrproc
!        j = rbuf(i+jproc)
!        if (j > 0) then
!          gmdistype(jproc) = j
!        endif
!      enddo
!    enddo
!    !
!    ! -- deallocate work arrays
!    deallocate(sbuf, rbuf, rcnt, offsets)
!    !
!    ! -- return
!    return
!  end subroutine mpi_get_dis_world

end module MpiExchangeModule
