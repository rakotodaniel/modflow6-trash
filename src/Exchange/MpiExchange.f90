module MpiExchangeModule
  use MpiExchangeGenModule, only: serialrun, writestd, partstr
  use KindModule, only: DP, I4B
  use SimModule, only: ustop, store_error
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, LENORIGIN,          &
                             LENVARNAME, LINELENGTH
  use gwfModule, only: gwfModelType
  use ListModule, only: ListType
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize, mpiwrpnrproc, mpiwrpmyrank,&
                        mpiwrpbarrier, mpiwrpcommworld, mpiwrpcomm_size,       &
                        mpiwrpcomm_rank, mpiwrpallgather, mpiwrpallgatherv,    &
                        mpiwrpcolstruct, mpiwrptypefree, ColMemoryType
  
  implicit none
  
  private
  
  ! -- Public types
  public :: MpiExchangeType
  ! -- Public variables
  public :: MpiWorld
  ! -- Public functions
  public :: mpi_initialize_world
  public :: mpi_dis_world
  public :: mpi_set_dis_world
    
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
    integer                                                :: lnmodel = 0   ! number of local models
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
    type(MpiGwfCommInt), dimension(:), pointer :: lxch      => null() ! MPI data structure point-to-point communication
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
    use MemoryManagerModule, only: mem_allocate !, mem_get_info
    ! -- dummy
    ! -- local
    character(len=LENORIGIN) :: origin
    integer(I4B) :: memitype, isize
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
    !serialrun = .false.
    ! @@@@@@ DEBUG
    !
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
    if (serialrun) then
      return
    endif
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
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(MpiExchangeType) :: this
    ! -- local
    character(len=LENORIGIN) :: origin
    class(NumericalModelType), pointer :: mp
    class(NumericalExchangeType), pointer :: cp
    integer(I4B) :: nm, im, ic, i, isub1, isub2, nja
    type(GwfModelType), pointer :: m1, m2
    character(len=LINELENGTH) :: errmsg
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    !
    ! -- Allocate scalars
    origin = this%name
    call mem_allocate(this%comm,   'COMM',   origin)
    call mem_allocate(this%nrproc, 'NRPROC', origin)
    call mem_allocate(this%myrank, 'MYRANK', origin)
    call mem_allocate(this%myproc, 'MYPROC', origin)
    call mem_allocate(this%nrxp,   'NRXP',   origin)
    !
    ! TODO: create communicator group; for now same as MPI_WORLD
    this%comm   = MpiWorld%comm
    this%nrproc = MpiWorld%nrproc
    this%myrank = MpiWorld%myrank
    this%myproc = MpiWorld%myproc
    !
    ! -- loop over my local models within this solution
    nm = this%lmodellist%Count()
    do im=1,nm
      mp => GetNumericalModelFromList(this%lmodellist, im)
      ! -- add model to local MPI model list
      call this%mpi_addmodel(2, mp%name)
      ! -- find the global subdomain number for this model
      i = ifind(MpiWorld%gmodelnames, mp%name)
      isub1 = MpiWorld%gsubs(i)
      ! -- add subdomain to local MPI subdomain list
      call this%mpi_addsub(2, isub1)
    enddo
    !
    ! -- loop over exchanges and count exchange partners within this solution
    this%nrxp = 0
    do ic=1,this%lexchangelist%Count()
      cp => GetNumericalExchangeFromList(this%lexchangelist, ic)
      call cp%get_m1m2(m1, m2)
      ! -- Add to interface
      if (cp%m2_ishalo) then
        this%nrxp = this%nrxp + 1
        ! -- Get my model subdomain number and check
        i = ifind(MpiWorld%gmodelnames, m1%name)
        isub1 = MpiWorld%gsubs(i)
        if (isub1 /= this%myproc) then
          write(errmsg,'(a)') 'Program error in mpi_local_exchange_init.'
          call store_error(errmsg)
          call ustop()
        endif
        i = ifind(MpiWorld%gmodelnames, m2%name)
        isub2 = MpiWorld%gsubs(i)
      end if
    end do
    !
    ! -- Allocate local communication data structure
    if (this%nrxp > 0) then
      allocate(this%lxch(this%nrxp))
    endif
    
    ! -- loop over exchanges and initialize
    this%nrxp = 0
    do ic=1,this%lexchangelist%Count()
      cp => GetNumericalExchangeFromList(this%lexchangelist, ic)
      call cp%get_m1m2(m1, m2)
      ! -- Add to interface
      if (cp%m2_ishalo) then
        this%nrxp = this%nrxp + 1
        i = ifind(MpiWorld%gmodelnames, m2%name)
        isub2 = MpiWorld%gsubs(i)
        ! -- allocate and set local exchange partner
        allocate(this%lxch(this%nrxp)%xprnk)
        this%lxch(this%nrxp)%xprnk = isub2 - 1
        ! --  allocate and set pointer to sending model (m1)
        allocate(this%lxch(this%nrxp)%sndgwfmodel1)
        this%lxch(this%nrxp)%sndgwfmodel1 => m1
        ! --  allocate and set pointer to receiving model (m2)
        allocate(this%lxch(this%nrxp)%recgwfmodel2) 
        this%lxch(this%nrxp)%recgwfmodel2 => m2
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

  subroutine mpi_dis_world(iopt)
! ******************************************************************************
! This subroutine gathers DIS information for all models.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeColModule, only: ciopt, n_recv, cmt_recv
    use MemoryManagerModule, only: mem_get_ptr, mem_setval
    use MemoryTypeModule, only: MemoryType
    use BaseModelModule, only: BaseModelType, GetBaseModelFromList  
    use ListsModule, only: basemodellist
    use NumericalModelModule, only: NumericalModelType
    ! -- dummy
    integer, intent(in) :: iopt
    ! -- local
    integer, parameter :: idis  = 1 
    integer, parameter :: idisv = 2 
    integer, parameter :: idisu = 3 
    !
    character(len=1) :: cdum
    character(len=4) :: dis_type
    integer :: im, is, iact
    class(BaseModelType), pointer :: mb 
    class(NumericalModelType), pointer :: mp
    type(MemoryType), pointer :: mt
    !
    integer, parameter :: cntypes = 6
    integer :: ierr, i, n_send, newtype
    integer, dimension(1) :: iwrk
    integer, dimension(:), allocatable :: recvcounts, displs
    type(ColMemoryType), dimension(:), allocatable :: cmt_send
    
    integer :: rank
! ------------------------------------------------------------------------------
    if (serialrun) then
      !return !@@@@@ DEBUG 
    endif
    
    do iact = 1, 2
      is = 0
      do im = 1, basemodellist%Count()
        mb => GetBaseModelFromList(basemodellist, im)
        select type (mb)
        class is (NumericalModelType)
          mp => mb
        end select  
        read(mp%dis%origin,*) cdum, dis_type
        if (iact == 1) then
          select case (dis_type)
          case ('DIS', 'DISV')
            if (iopt == 1) then
              is = is + 1
            else  
              is = is + 4
            endif
          case ('DISU')
            if (iopt == 1) then
              is = is + 1
            else
              is = is + 2
            endif
          end select
        else
        select case (dis_type)
          case ('DIS')
            if (iopt == 1) then
              !call mem_get_ptr('NEQ', 'GWF_MODEL_1', mt)
              !mp%neq = 456
              !call mem_get_ptr('NEQ', 'GWF_MODEL_1', mt)
              !call mem_setval(789,'NEQ','GWF_MODEL_1')
              !call mem_get_ptr('NEQ', 'GWF_MODEL_1', mt)
              !mp%neq = 456
              !call mem_get_ptr('NEQ', 'GWF_MODEL_1', mt)
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            else
              call mem_get_ptr('NLAY', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('NROW', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('NCOL', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            endif
          case ('DISV')
            if (iopt == 1) then
              call mem_get_ptr('DNDIM',  mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            else
              call mem_get_ptr('NLAY',  mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('NCPL',  mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('NVERT', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            endif
          case ('DISU')
            if (iopt == 1) then
              call mem_get_ptr('DNDIM', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            else  
              call mem_get_ptr('NVERT', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
              call mem_get_ptr('MSHAPE', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            endif
          end select
        endif
      enddo
      if (iact == 1) then
        n_send = is
        if (n_send > 0) then
          allocate(cmt_send(n_send))
        endif
      endif
    enddo
    !
    ! -- gather sizes
    allocate(recvcounts(MpiWorld%nrproc), displs(MpiWorld%nrproc))
    iwrk(1) = n_send
    call mpiwrpallgather(MpiWorld%comm, iwrk, 1, recvcounts, 1)
    n_recv = sum(recvcounts)
    if (allocated(cmt_recv)) then
      deallocate(cmt_recv)
    endif
    allocate(cmt_recv(n_recv))
    displs(1) = 0
    do i = 2, MpiWorld%nrproc
      displs(i) = displs(i-1) + recvcounts(i-1)
    end do
    ! -- create derived type
    call mpiwrpcolstruct(newtype)
    ! -- gather data
    call mpiwrpallgatherv(MpiWorld%comm, cmt_send, n_send, newtype, cmt_recv,   &
                          recvcounts, newtype, displs)
    ! -- DEBUG
    do rank = 0, MpiWorld%nrproc-1
      if (MpiWorld%myrank == rank) then
        write(*,*) '=== send myrank', rank
        do i = 1, n_send
          write(*,'(a,1x,a,1x,i)') trim(cmt_send(i)%name), trim(cmt_send(i)%origin), cmt_send(i)%intsclr  
        end do
        write(*,*) '=== recv myrank', rank
        do i = 1, n_recv
          write(*,'(a,1x,a,1x,i)') trim(cmt_recv(i)%name), trim(cmt_recv(i)%origin), cmt_recv(i)%intsclr 
        end do
      endif
      call mpiwrpbarrier(MpiWorld%comm)
    enddo
    !
    ! -- clean up
    call mpiwrptypefree(newtype)
    deallocate(cmt_send, recvcounts, displs)
    ciopt = iopt
    !
    ! -- return
    return
  end subroutine mpi_dis_world

  subroutine mpi_set_dis_world()
! ******************************************************************************
! This subroutine sets the DIS scalars for the halo (m2) models.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeColModule, only: ciopt, n_recv, cmt_recv
    use MemoryTypeModule, only: ilogicalsclr, iintsclr, idblsclr,               & 
                                iaint1d, iaint2d,                               & 
                                iadbl1d, iadbl2d
    use MemoryManagerModule, only: mem_setval
    use ArrayHandlersModule, only: ExpandArray, ifind
    use BaseExchangeModule, only: BaseExchangeType, GetBaseExchangeFromList
    use NumericalExchangeModule, only: NumericalExchangeType
    use ListsModule, only: baseexchangelist
    ! -- dummy
    ! -- local
    class(BaseExchangeType),   pointer :: bep
    class(NumericalExchangeType), pointer :: nep
    type(GwfModelType), pointer :: m1, m2
    character(len=LENVARNAME)   :: name        !name of the array
    character(len=LENORIGIN)    :: origin      !name of origin
    character(len=LENMODELNAME) :: modelname   !name of origin
    integer :: i, ic, n, m
    character(len=LENMODELNAME), allocatable, dimension(:) :: modelname_halo
! ------------------------------------------------------------------------------
    !
    n = 0
    do ic=1,baseexchangelist%Count()
      bep => GetBaseExchangeFromList(baseexchangelist, ic)
      select type (bep)
      class is (NumericalExchangeType)
        nep => bep
      end select      
      call nep%get_m1m2(m1, m2)
      if (nep%m2_ishalo) then
        m = ifind(modelname_halo, m2%name)
        if (m < 0) then
          n = n + 1
          call ExpandArray(modelname_halo)
          modelname_halo(n) = m2%name
        endif
      endif
    enddo 
    !
    do i = 1, n_recv
      name   = cmt_recv(i)%name
      origin = cmt_recv(i)%origin
      read(origin,*) modelname
      m = ifind(modelname_halo, modelname)
      if (m > 0) then
        select case(cmt_recv(i)%memitype)
          case(iintsclr)
            write(*,*) 'Setting integer for model '//trim(modelname)//' '//trim(name)//' '//trim(origin)
            call mem_setval(cmt_recv(i)%intsclr, name, origin)
          case(iaint1d)
            if (MpiWorld%myrank == 0) then
              write(*,*) MpiWorld%myrank
              write(*,*) 'Setting array for model '//trim(modelname)//' '//trim(name)//' '//trim(origin)
            endif
            call mem_setval(cmt_recv(i)%aint1d, name, origin)
        end select
      endif
    end do
    !
    ! -- cleanup
    deallocate(modelname_halo)
    deallocate(cmt_recv)
  
    ! -- return
    return
  end subroutine mpi_set_dis_world
  
  subroutine mpi_to_colmem(mt, cmt)
! ******************************************************************************
! Convert to collective MemoryType
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use ConstantsModule, only: DZERO  
    use MemoryTypeModule, only: MemoryType
    ! -- dummy
    type(MemoryType), intent(in) :: mt
    type(ColMemoryType), intent(out) :: cmt
    ! -- local
! ------------------------------------------------------------------------------
    cmt%name     = mt%name
    cmt%origin   = mt%origin
    cmt%memitype = mt%memitype
    if (associated(mt%logicalsclr)) then
      cmt%logicalsclr = mt%logicalsclr
    else
      cmt%logicalsclr = .false.
    endif  
    if (associated(mt%intsclr)) then
      cmt%intsclr = mt%intsclr
    else
      cmt%intsclr = 0
    endif
    if (associated(mt%dblsclr)) then
      cmt%dblsclr = mt%dblsclr
    else
      cmt%dblsclr = DZERO
    endif
    if (associated(mt%aint1d)) then
      cmt%aint1d = mt%aint1d
    else
      cmt%aint1d = 0
    endif
    ! -- return
    return
  end subroutine mpi_to_colmem  
end module MpiExchangeModule
