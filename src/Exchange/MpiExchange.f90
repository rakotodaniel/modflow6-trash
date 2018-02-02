module MpiExchangeModule
  
  use KindModule, only: DP, I4B  
  use SimModule, only: ustop, store_error
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, LENORIGIN,          &
                             LENVARNAME, LINELENGTH
  use ArrayHandlersModule, only: ExpandArray, ifind
  use gwfModule, only: gwfModelType
  use ListModule, only: ListType
  use MpiExchangeGenModule, only: serialrun, writestd, partstr
  use MpiWrapper, only: mpiwrpinit, mpiwrpfinalize, mpiwrpnrproc, mpiwrpmyrank,&
                        mpiwrpbarrier, mpiwrpcommworld, mpiwrpcomm_size,       &
                        mpiwrpcomm_rank, mpiwrpallgather, mpiwrpallgatherv,    &
                        mpiwrpcolstruct, mpiwrptypefree, ColMemoryType,        &
                        mpiwrpstats, mpiwrpisend, mpiwrpirecv,                 &
                        mpiwrpmmtstruct, mpiwrpmtstruct, mpiwrpwaitall,        &
                        mpiwrpprobe, mpiwrpgetcount
  use MemoryTypeModule, only: MemoryType
  use MpiWrapper, only: MetaMemoryType
  use NumericalExchangeModule, only: NumericalExchangeType

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
  public :: mpi_add_halo_model

  ! -- Container for NumericalExchangeType
  type NumericalExchangeContType
    class(NumericalExchangeType), pointer :: obj => null()
  end type NumericalExchangeContType
    
  ! -- Local types
  type :: MpiGwfCommInt
    integer(I4B), pointer                                  :: xprnk        => null() ! neighboring rank
    integer(I4B), pointer                                  :: nexchange    => null() ! number of exchanges 
    type(NumericalExchangeContType), dimension(:), pointer :: exchange     => null() ! exchanges
    character(len=LENPACKAGENAME), dimension(:), pointer   :: name         => null() ! variable name
    character(len=LENORIGIN), dimension(:), pointer        :: origin       => null() ! variable origin
    integer(I4B)                                           :: nsndvar = 0            ! number of variables to send
    integer(I4B)                                           :: nrcvvar = 0            ! number of variables to receive
    type(MemoryType), dimension(:), allocatable            :: sndmt                  ! send data
    type(MemoryType), dimension(:), allocatable            :: rcvmt                  ! receive data
    type(MetaMemoryType), dimension(:), allocatable        :: sndmmt                 ! send meta data
    type(MetaMemoryType), dimension(:), allocatable        :: rcvmmt                 ! receive meta data
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
    procedure :: mpi_local_exchange
    procedure :: mpi_debug
    procedure :: mpi_da
  end type MpiExchangeType
  
  ! -- World communicator
  type(MpiExchangeType), pointer :: MpiWorld => null()

  character(len=LINELENGTH) :: errmsg

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
    character(len=LENMODELNAME) :: mname
    class(NumericalModelType), pointer :: mp
    class(NumericalExchangeType), pointer :: cp
    integer(I4B) :: nm, im, ic, i, ip, ixp, isub1, isub2, nja, nex
    type(GwfModelType), pointer :: m1, m2
    integer(I4B), dimension(:), allocatable :: iwrk
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return !@@@@@DEBUG
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
      if (mp%ishalo) then
        read(mp%name,*) mname
      else
        mname = mp%name
      endif
      ! -- find the global subdomain number for this model
      i = ifind(MpiWorld%gmodelnames, mname)
      ! -- add subdomain to local MPI subdomain list
      if (i > 0) then
        isub1 = MpiWorld%gsubs(i)
        call this%mpi_addsub(2, isub1)
      else
        write(errmsg,'(a)') 'Program error in mpi_local_exchange_init.'
        call store_error(errmsg)
        call ustop()
      endif
    enddo
    !
    ! -- loop over exchanges and count exchange partners within this solution
    allocate(iwrk(this%nrproc))
    iwrk = 0
    do ic=1,this%lexchangelist%Count()
      cp => GetNumericalExchangeFromList(this%lexchangelist, ic)
      call cp%get_m1m2(m1, m2)
      ! -- Add to interface
      if (cp%m2_ishalo) then
        ! -- Get my model subdomain number and check
        i = ifind(MpiWorld%gmodelnames, m1%name)
        isub1 = MpiWorld%gsubs(i)
        if (isub1 /= this%myproc) then
          write(errmsg,'(a)') 'Program error in mpi_local_exchange_init.'
          call store_error(errmsg)
          call ustop()
        endif
        read(m2%name,*) mname
        i = ifind(MpiWorld%gmodelnames, mname)
        isub2 = MpiWorld%gsubs(i)
        iwrk(isub2) = iwrk(isub2) + 1
      end if
    end do
    this%nrxp = 0
    do ip = 1, this%nrproc
      if (iwrk(ip) > 0) then
        this%nrxp = this%nrxp + 1
      endif
    enddo
    !
    ! -- Allocate local communication data structure
    if (this%nrxp > 0) then
      allocate(this%lxch(this%nrxp))
    endif
    !
    ixp = 0
    do ip = 1, this%nrproc
      nex = iwrk(ip)
      if (nex > 0) then
        ixp = ixp + 1
        allocate(this%lxch(ixp)%nexchange)
        allocate(this%lxch(ixp)%exchange(nex))
        allocate(this%lxch(ixp)%xprnk)
        this%lxch(ixp)%nexchange = 0
        this%lxch(ixp)%xprnk = ip-1
      endif
      ! -- Set mapping to exchange partner index
      iwrk(ip) = ixp
    enddo
    !
    ! -- loop over exchanges and initialize
    do ic=1,this%lexchangelist%Count()
      cp => GetNumericalExchangeFromList(this%lexchangelist, ic)
      call cp%get_m1m2(m1, m2)
      ! -- Add to interface
      if (cp%m2_ishalo) then
        ! -- Get my model subdomain number and check
        read(m2%name,*) mname
        i = ifind(MpiWorld%gmodelnames, mname)
        isub2 = MpiWorld%gsubs(i)
        ixp = iwrk(isub2)
        nex = this%lxch(ixp)%nexchange
        nex = nex + 1
        ! -- set pointer to exchange
        this%lxch(ixp)%exchange(nex)%obj => cp
        this%lxch(ixp)%nexchange = nex
      end if
    end do
    
    ! -- Clean up
    deallocate(iwrk)
    !
    ! -- return
    return
  end subroutine mpi_local_exchange_init

  subroutine mpi_local_exchange(this, solname, iopt)
! ******************************************************************************
! Local point-to-point exchange (wrapper).
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use NumericalExchangeModule, only: NumericalExchangeType,                  &
                                    GetNumericalExchangeFromList
    use GwfModule, only: GwfModelType
    use MemoryManagerModule, only: mem_get_ptr, mem_setval
    use MpiWrapper, only: mpiwrpstats
    ! -- dummy
    class(MpiExchangeType) :: this
    character(len=*) :: solname
    integer(I4B), intent(in) :: iopt
    ! -- local
    integer, parameter :: MAXNVAR = 100
    type VarMemoryType
      integer                                       :: nvar = 0
      character(len=LENVARNAME), dimension(MAXNVAR) :: names
      character(len=LENVARNAME), dimension(MAXNVAR) :: namesext
      logical, dimension(MAXNVAR)                   :: namessol
    end type VarMemoryType
    type(VarMemoryType), pointer :: vmt1, vmt2, vmt
    
    integer(I4B), parameter :: nvar1 = 2
    character(len=LENVARNAME), dimension(nvar1) :: names1, namesext1
    logical, dimension(nvar1) :: namessol1
    data names1    /'AREA',  'IACTIVE'/
    data namesext1 /'DIS',   ''/
    data namessol1 /.false., .true./
    !
    class(NumericalExchangeType), pointer :: cp
    type(GwfModelType), pointer :: m1, m2
    type(MemoryType), pointer :: mt
    character(len=LENORIGIN) :: mod_origin, sol_origin, src_origin, tgt_origin
    integer(I4B) :: ixp, ix, nsnd, iv, is, i, istat
    character(len=LENMODELNAME) :: mname, id
    !
    integer(I4B) :: newtype
    integer(I4B), dimension(:), allocatable :: snd_newtype, rcv_newtype
    integer(I4B), dimension(:), allocatable :: sreq, rreq
    !
    integer(I4B) :: irank !@@@DEBUG
! ------------------------------------------------------------------------------
    !
    if (serialrun) then
      return
    endif
    !
    ! -- Set the variabels
    allocate(vmt1)
    iv = 0
    iv = iv + 1
    vmt1%names(iv)    = 'TOP'
    vmt1%namesext(iv) = 'DIS'
    vmt1%namessol(iv) = .false.
    iv = iv + 1
    vmt1%names(iv)    = 'BOT'
    vmt1%namesext(iv) = 'DIS'
    vmt1%namessol(iv) = .false.
    iv = iv + 1
    vmt1%names(iv)    = 'AREA'
    vmt1%namesext(iv) = 'DIS'
    vmt1%namessol(iv) = .false.
    iv = iv + 1
    vmt1%names(iv)    = 'IACTIVE'
    vmt1%namesext(iv) = ''
    vmt1%namessol(iv) = .true.
    vmt1%nvar = iv
    !
    if (iopt == 1) then
      vmt => vmt1
    endif
    !
    ! -- Loop over exchange partners
    do ixp = 1, this%nrxp
      nsnd = this%lxch(ixp)%nexchange * vmt%nvar
      this%lxch(ixp)%nsndvar = nsnd
      ! -- Allocate send buffers
      allocate(this%lxch(ixp)%sndmt(nsnd))
      allocate(this%lxch(ixp)%sndmmt(nsnd))
      is = 0
      do ix = 1, this%lxch(ixp)%nexchange
        cp => this%lxch(ixp)%exchange(ix)%obj 
        call cp%get_m1m2(m1, m2)
        do iv = 1, vmt%nvar
          write(mod_origin,'(a,1x,a)') trim(m1%name), trim(vmt%namesext(iv))
          sol_origin = trim(solname)
          if (.not.vmt%namessol(iv)) then
            src_origin = mod_origin
          else
            src_origin = sol_origin
          endif
          tgt_origin = mod_origin
          call mem_get_ptr(vmt%names(iv), src_origin, mt)
          is = is + 1
          call mpi_pack_mt(tgt_origin, cp%nodem1, cp%nexg, mt, this%lxch(ixp)%sndmt(is))
          call mpi_pack_mmt(this%lxch(ixp)%sndmt(is), this%lxch(ixp)%sndmmt(is))
        enddo
      enddo
    enddo
    !    
    ! -- First point-to-point communication to exchange memory type meta-data
    !
    ! -- Create MPI derived datatype
    call mpiwrpmmtstruct(newtype)
    ! -- Allocate request arrays
    allocate(sreq(this%nrproc), rreq(this%nrproc))
    ! -- Non-blocking send
    do ixp = 1, this%nrxp
      call mpiwrpisend(this%lxch(ixp)%sndmmt, this%lxch(ixp)%nsndvar, newtype,  &
                       this%lxch(ixp)%xprnk, 0, this%comm, sreq(ixp))
    enddo
    call mpiwrpwaitall(this%nrproc, sreq, mpiwrpstats)
    ! -- Probe for the data sizes and allocate
    do ixp = 1, this%nrxp
       call mpiwrpprobe(this%lxch(ixp)%xprnk, 0, this%comm, mpiwrpstats(1,ixp))
       call mpiwrpgetcount(mpiwrpstats(1,ixp), newtype, this%lxch(ixp)%nrcvvar)
       allocate(this%lxch(ixp)%rcvmmt(this%lxch(ixp)%nrcvvar))
       allocate(this%lxch(ixp)%rcvmt (this%lxch(ixp)%nrcvvar))
    enddo
    ! -- Non-blocking receive
    do ixp = 1, this%nrxp
      call mpiwrpirecv(this%lxch(ixp)%rcvmmt, this%lxch(ixp)%nrcvvar, newtype,  &
                       this%lxch(ixp)%xprnk, 0, this%comm, rreq(ixp))
    enddo
    call mpiwrpwaitall(this%nrproc, rreq, mpiwrpstats)
    ! -- Clean up MPI datatype
    call mpiwrptypefree(newtype)
    ! -- Debug
    if (.false.) then
      do irank = 0, this%nrproc-1
        if (irank == this%myrank) then
          write(*,*) '====myrank',this%myrank
          do ixp = 1, this%nrxp
            write(*,*) '=======received from',this%lxch(ixp)%xprnk
            do i = 1, this%lxch(ixp)%nrcvvar
              write(*,*) 'name:     ', trim(this%lxch(ixp)%rcvmmt(i)%name)
              write(*,*) 'origin:   ', trim(this%lxch(ixp)%rcvmmt(i)%origin)
              write(*,*) 'memitype: ', this%lxch(ixp)%rcvmmt(i)%memitype
              write(*,*) 'isize:    ', this%lxch(ixp)%rcvmmt(i)%isize
              write(*,*) 'ncol:     ', this%lxch(ixp)%rcvmmt(i)%ncol
              write(*,*) 'nrow:     ', this%lxch(ixp)%rcvmmt(i)%nrow
            enddo
          enddo
        endif
        call mpiwrpbarrier(this%comm)
      enddo
    endif
    !
    ! -- Second point-to-point communication to exchange actual data
    !
    ! -- Create MPI send and receive datatypes
    allocate(snd_newtype(this%nrxp), rcv_newtype(this%nrxp))
    do ixp = 1, this%nrxp
      call mpiwrpmtstruct(this%lxch(ixp)%sndmt, this%lxch(ixp)%sndmmt,          &
                          this%lxch(ixp)%nsndvar, snd_newtype(ixp))
      call mpiwrpmtstruct(this%lxch(ixp)%rcvmt, this%lxch(ixp)%rcvmmt,          &
                          this%lxch(ixp)%nrcvvar, rcv_newtype(ixp))
    enddo
    ! -- Non-blocking send
    do ixp = 1, this%nrxp
      call mpiwrpisend(this%lxch(ixp)%sndmt, 1, snd_newtype(ixp),               &
                       this%lxch(ixp)%xprnk, 0, this%comm, sreq(ixp))
    enddo
    call mpiwrpwaitall(this%nrproc, sreq, mpiwrpstats)
    ! -- Non-blocking receive
    do ixp = 1, this%nrxp
      call mpiwrpirecv(this%lxch(ixp)%rcvmt, 1, rcv_newtype(ixp),               &
                       this%lxch(ixp)%xprnk, 0, this%comm, rreq(ixp))
    enddo
    call mpiwrpwaitall(this%nrproc, rreq, mpiwrpstats)
    ! -- Clean up MPI datatypes
    do ixp = 1, this%nrxp
      call mpiwrptypefree(snd_newtype(ixp))
      call mpiwrptypefree(rcv_newtype(ixp))
    enddo
    deallocate(snd_newtype, rcv_newtype)
    ! -- Debug
    if (.false.) then
      do irank = 0, this%nrproc-1
        if (irank == this%myrank) then
          write(*,*) '====myrank',this%myrank
          do ixp = 1, this%nrxp
            write(*,*) '=======received from',this%lxch(ixp)%xprnk
            do i = 1, this%lxch(ixp)%nrcvvar
              write(*,*) 'name:      ', trim(this%lxch(ixp)%rcvmt(i)%name)
              write(*,*) 'origin:    ', trim(this%lxch(ixp)%rcvmt(i)%origin)
              write(*,*) 'isize:     ', this%lxch(ixp)%rcvmt(i)%isize
              write(*,*) 'memitype:  ', this%lxch(ixp)%rcvmt(i)%memitype
              !write(*,*) 'dval:      ', this%lxch(ixp)%rcvmt(i)%adbl1d(1)
            enddo
          enddo
        endif
        call mpiwrpbarrier(this%comm)
      enddo
    endif
    ! 
    ! -- Set the received data for the halo (m2) models
    do ixp = 1, this%nrxp
      do ix = 1, this%lxch(ixp)%nexchange
        cp => this%lxch(ixp)%exchange(ix)%obj 
        call cp%get_m1m2(m1, m2)
        do i = 1, this%lxch(ixp)%nrcvvar
          src_origin = this%lxch(ixp)%rcvmt(i)%origin
          read(src_origin,*,iostat=istat) mname, id
          if (istat /= 0) then
            read(src_origin,*,iostat=istat) mname
            id = ''
          endif
          write(tgt_origin, '(a,1x,a)') trim(m2%name), trim(id)
          this%lxch(ixp)%rcvmt(i)%origin = tgt_origin
          call mem_setval(this%lxch(ixp)%rcvmt(i))
        enddo
      enddo 
    enddo
    !
    ! -- Debug
    !call ustop(stopmess='@@@@ STOP LOCAL EXCHANGE')
    !
    ! -- return
    return
  end subroutine mpi_local_exchange

  subroutine mpi_pack_mt(origin, node, nexg, mti, mto)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iaint1d, iadbl1d
    ! -- dummy
    character(len=*), intent(in) :: origin
    integer(I4B), intent(in) :: nexg
    integer(I4B), dimension(nexg), intent(in) :: node
    type(MemoryType), intent(in) :: mti
    type(MemoryType), intent(out) :: mto
    ! -- local
    integer(I4B) :: i, n
! ------------------------------------------------------------------------------
    !
    write(errmsg,'(a)') 'Program error in mpi_pack_mt.'
    !
    mto%name     = mti%name
    mto%origin   = trim(origin)
    mto%memitype = mti%memitype
    mto%isize    = nexg
    !
    select case(mti%memitype)
      case(iaint1d)
        allocate(mto%aint1d(nexg))
        do i = 1, nexg
          n = node(i)
          if (n < 0 .or. n > size(mti%aint1d)) then
            call store_error(errmsg)
            call ustop()
          endif
          mto%aint1d(i) = mti%aint1d(n)
        enddo
        !write(*,*) '# int n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%aint1d)
      case(iadbl1d)
        allocate(mto%adbl1d(nexg))
        do i = 1, nexg
          n = node(i)
          if (n < 0 .or. n > size(mti%adbl1d)) then
            call store_error(errmsg)
            call ustop()
          endif
          mto%adbl1d(i) = mti%adbl1d(n)
        enddo
        !write(*,*) '# dbl n, isize',trim(mti%name)//' '//trim(mti%origin), n,size(mti%adbl1d)
      case default
        call store_error(errmsg)
        call ustop()
    end select
    !
    ! -- return
    return
  end subroutine mpi_pack_mt
  
  subroutine mpi_pack_mmt(mt, mmt)
! ******************************************************************************
! Pack memory type for point-to-point communication.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MemoryTypeModule, only: iaint2d, iadbl2d
    ! -- dummy
    type(MemoryType), intent(in) :: mt
    type(MetaMemoryType), intent(out) :: mmt
    ! -- local
! ------------------------------------------------------------------------------
    !
    mmt%name     = mt%name
    mmt%origin   = mt%origin
    mmt%memitype = mt%memitype
    mmt%isize    = mt%isize
    select case(mt%memitype)
    case(iaint2d)
        mmt%ncol = size(mt%aint2d, dim=1)
        mmt%nrow = size(mt%aint2d, dim=2)
      case(iadbl2d)
        mmt%ncol = size(mt%adbl2d, dim=1)
        mmt%nrow = size(mt%adbl2d, dim=2)
    end select
    !
    ! -- return
    return
  end subroutine mpi_pack_mmt

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
    integer(I4B), intent(in)  :: isub
    logical, intent(out) :: add
    ! -- local
! ------------------------------------------------------------------------------
    if (serialrun) then
      !add = .true. !DEBUG 
      !return
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
! ------------------------------------------------------------------------------
    if (serialrun) then
      ! return @@@ DEBUG
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
              is = is + 7
            endif
          case ('DISU')
            if (iopt == 1) then
              is = is + 1
            else
              is = is + 5
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
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
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
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
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
              call mem_get_ptr('NODESUSER', mp%dis%origin, mt)
              is = is + 1
              call mpi_to_colmem(mt, cmt_send(is))
            endif
          end select
          if (iopt == 2) then
            call mem_get_ptr('INNPF', mp%name, mt)
            is = is + 1
            call mpi_to_colmem(mt, cmt_send(is))
            call mem_get_ptr('INIC', mp%name, mt)
            is = is + 1
            call mpi_to_colmem(mt, cmt_send(is))
          endif
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
    if (.false.) then
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
    endif
    !
    ! -- clean up
    call mpiwrptypefree(newtype)
    deallocate(cmt_send, recvcounts, displs)
    ciopt = iopt
    !
    ! -- return
    return
  end subroutine mpi_dis_world
    
  subroutine mpi_add_halo_model(im, modelname)
! ******************************************************************************
! This subroutine sets the list of halo (m2) models
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    use MpiExchangeGenModule, only: nhalo, modelname_halo,                      &
                                    mpi_create_modelname_halo
    ! -- dummy
    integer, intent(in) :: im
    character(len=*), intent(inout) :: modelname
    ! -- local
    integer(I4B) :: m
! ------------------------------------------------------------------------------
    call mpi_create_modelname_halo(im, modelname)
    m = ifind(modelname_halo, modelname)
    if (m < 0) then
      nhalo = nhalo + 1
      call ExpandArray(modelname_halo)
      modelname_halo(nhalo) = modelname
    endif
    !
    ! -- return
    return
  end subroutine mpi_add_halo_model
  
  subroutine mpi_set_dis_world()
! ******************************************************************************
! This subroutine sets the DIS scalars for the halo (m2) models.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    ! -- modules
    use MpiExchangeGenModule, only: modelname_halo, nhalo,                      &
                                    mpi_create_modelname_halo
    use MpiExchangeColModule, only: ciopt, n_recv, cmt_recv
    use MemoryTypeModule, only: ilogicalsclr, iintsclr, idblsclr,               & 
                                iaint1d, iaint2d,                               & 
                                iadbl1d, iadbl2d
    use MemoryManagerModule, only: mem_setval
    ! -- dummy
    ! -- local
    character(len=LENVARNAME)   :: name        !name of the array
    character(len=LENORIGIN)    :: origin      !name of origin
    character(len=LENORIGIN)    :: origin_halo !name of origin
    character(len=LENMODELNAME) :: mname       !name of origin
    character(len=LENMODELNAME) :: mname_halo  !name of origin
    character(len=LENMODELNAME) :: id
    integer(I4B) :: i, j, m, istat
    
    integer(I4B) :: wrk(3) !DEBUG
! ------------------------------------------------------------------------------
    !
    ! -- Set the list of halo model names
    !call mpi_set_modelname_halo()
    !
    if (serialrun) then
      !call mem_setval(1, 'NLAY', 'GWF_MODEL_2 HALO1 DIS') !@@@DEBUG
      !call mem_setval(3, 'NROW', 'GWF_MODEL_2 HALO1 DIS')
      !call mem_setval(3, 'NCOL', 'GWF_MODEL_2 HALO1 DIS')
      !wrk(1) = 1
      !wrk(2) = 3
      !wrk(3) = 3
      !call mem_setval(wrk, 'MSHAPE', 'GWF_MODEL_2 HALO1 DIS')
      !call mem_setval(9, 'NODESUSER', 'GWF_MODEL_2 HALO1 DIS')
      !call mem_setval(1, 'INNPF', 'GWF_MODEL_2 HALO1')
      !call mem_setval(1, 'INIC', 'GWF_MODEL_2 HALO1')
      return
    endif
    !
    do i = 1, n_recv
      name   = cmt_recv(i)%name
      origin = cmt_recv(i)%origin
      read(origin,*,iostat=istat) mname, id
      if (istat /= 0) then
        read(origin,*,iostat=istat) mname
        id = ''
      endif
      do j = 1, nhalo
        mname_halo = mname
        call mpi_create_modelname_halo(j, mname_halo)
        m = ifind(modelname_halo, mname_halo)
        if (m > 0) then
          origin_halo = trim(mname_halo)//' '//trim(id)
          select case(cmt_recv(i)%memitype)
          case(iintsclr)
              if (MpiWorld%myrank == 1) then
                !write(*,*) MpiWorld%myrank
                !write(*,*) 'Setting integer: '//trim(name)//', "'//trim(origin_halo)//'"', cmt_recv(i)%intsclr
              endif
              call mem_setval(cmt_recv(i)%intsclr, name, origin_halo)
            case(iaint1d)
              if (MpiWorld%myrank == 1) then
                !write(*,*) MpiWorld%myrank
                !write(*,*) 'Setting array for model '//trim(mname_halo)//': '//trim(name)//', '//trim(origin_halo)
              endif
              call mem_setval(cmt_recv(i)%aint1d, name, origin_halo)
          end select
        endif
      enddo
    end do
    !
    ! -- cleanup
    deallocate(cmt_recv)
    !
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