!
program mf6
! ******************************************************************************
! Main MODFLOW Version 6 program.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
  ! -- modules
  use KindModule,             only: DP, I4B
  use ConstantsModule,        only: MFVNAM, VERSION, MFTITLE, FMTDISCLAIMER,   &
                                    ISTDOUT, IDEVELOPMODE
  use CompilerVersion
  use InputOutputModule,      only: write_centered
  use SimulationCreateModule, only: simulation_cr, simulation_da
  use TimerModule,            only: start_time, elapsed_time
  use MemoryManagerModule,    only: mem_usage, mem_da
  use BaseModelModule,        only: BaseModelType, GetBaseModelFromList
  use BaseExchangeModule,     only: BaseExchangeType, GetBaseExchangeFromList
  use BaseSolutionModule,     only: BaseSolutionType, GetBaseSolutionFromList
  use SolutionGroupModule,    only: SolutionGroupType, GetSolutionGroupFromList
  use ListsModule,            only: basesolutionlist, solutiongrouplist,       &
                                    basemodellist, baseexchangelist,           &
                                    lists_da,                                  &
                                    halomodellist !JV
  use SimVariablesModule,     only: iout 
  use SimModule,              only: converge_reset, converge_check,            &
                                    final_message
  use TdisModule,             only: tdis_tu, tdis_da,                          &
                                    endofsimulation
  use MpiExchangeGenModule,   only: mpi_initialize, serialrun, writestd !JV
  use MpiExchangeModule,      only: mpi_initialize_world, mpi_dis_world,       & !JV
                                    mpi_set_dis_world, MpiWorld
  use NumericalSolutionModule, only: NumericalSolutionType !JV
  implicit none
  ! -- local
  class(SolutionGroupType), pointer :: sgp
  class(BaseSolutionType),  pointer :: sp
  class(NumericalSolutionType), pointer :: nsp !JV
  class(BaseModelType),     pointer :: mp
  class(BaseExchangeType),  pointer :: ep
  integer(I4B) :: im, ic, is, isg
  logical :: exit_tsloop
  character(len=80) :: compiler
  ! -- formats
! ------------------------------------------------------------------------------
  !
  ! -- Initialize MPI if required
  call mpi_initialize() !JV
  call mpi_initialize_world() !JV
  !
  ! -- Write banner to screen (unit 6) and start timer
  call write_centered('MODFLOW'//MFVNAM, ISTDOUT, 80)
  call write_centered(MFTITLE, ISTDOUT, 80)
  call write_centered('VERSION '//VERSION, ISTDOUT, 80)
  if (.not.serialrun) call write_centered('***RUNNING IN PARALLEL MODE WITH '& !JV
    //TRIM(MpiWorld%nrprocstr)//' MPI PROCESSES***',ISTDOUT, 80) !JV
  !
  ! -- Write if develop mode
  if (IDEVELOPMODE == 1) call write_centered('***DEVELOP MODE***', ISTDOUT, 80)
  !
  ! -- Write compiler version
  call get_compiler(compiler)
  call write_centered(' ', ISTDOUT, 80)
  call write_centered(trim(adjustl(compiler)), ISTDOUT, 80)
  !
  ! -- Write disclaimer
  if (writestd) write(ISTDOUT, FMTDISCLAIMER) !JV
  ! -- get start time
  call start_time()
  !
  !
  ! -- CREATE (CR)
  call simulation_cr()
  !
  !
  ! -- DEFINE (DF)
  ! -- Define each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_df()
  enddo
  !
  ! -- Collective MPI communication scalars DIS
  call mpi_dis_world(2) !JV
  call mpi_set_dis_world() !JV
  !
  ! -- Define each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_df()
  enddo
  !
  ! -- Define each solution
  do is = 1, basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_df()
  enddo
  !
  !
  ! -- ALLOCATE AND READ (AR)
  ! -- Allocate and read each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_ar()
  enddo
  ! -- Allocate and read each model
  do im = 1, halomodellist%Count() !JV
    mp => GetBaseModelFromList(halomodellist, im) !JV
    call mp%model_ar() !JV
  enddo !JV
  ! 
  ! -- Local exchange (TODO)
  do is=1,basesolutionlist%Count() !JV
    sp => GetBaseSolutionFromList(basesolutionlist, is) !JV
    select type (sp) !JV
    class is (NumericalSolutionType) !JV
      nsp => sp !JV
    end select !JV
    call nsp%MpiSol%mpi_local_exchange(nsp%name, 1) !JV
  enddo !JV
  !
  ! -- Allocate and read each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_ar()
  enddo
  !
  ! -- Allocate and read each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_ar()
  enddo
  !
  !
  ! -- TIME STEP LOOP
  tsloop: do
    !
    ! -- TIME UPDATE (TU)
    call tdis_tu()
    !
    !
    ! -- READ AND PREPARE (RP)
    ! -- Read and prepare each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_rp()
    enddo
    !
    ! -- Read and prepare each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_rp()
    enddo
    !
    ! -- Read and prepare each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_rp()
    enddo
    !
    !
    ! -- CALCULATE (CA)
    call converge_reset()
    do isg = 1, solutiongrouplist%Count()
      sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
      call sgp%sgp_ca()
    enddo
    !
    !
    ! -- OUTPUT (OT)
    ! -- Write output for each model
    do im = 1, basemodellist%Count()
      mp => GetBaseModelFromList(basemodellist, im)
      call mp%model_ot()
    enddo
    !
    ! -- Write output for each exchange
    do ic = 1, baseexchangelist%Count()
      ep => GetBaseExchangeFromList(baseexchangelist, ic)
      call ep%exg_ot()
    enddo
    !
    ! -- Write output for each solution
    do is=1,basesolutionlist%Count()
      sp => GetBaseSolutionFromList(basesolutionlist, is)
      call sp%sln_ot()
    enddo
    !
    ! -- Time step exit conditions
    call converge_check(exit_tsloop)
    if(exit_tsloop) exit tsloop
    if(endofsimulation) exit tsloop
    !
  enddo tsloop
  !
  !
  ! -- FINAL PROCESSING (FP)
  ! -- Final processing for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_fp()
  enddo
  !
  ! -- Final processing for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_fp()
  enddo
  !
  ! -- Final processing for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_fp()
  enddo
  !
  ! -- DEALLOCATE (DA)
  ! -- Deallocate tdis
  call tdis_da()
  !
  ! -- Deallocate for each model
  do im = 1, basemodellist%Count()
    mp => GetBaseModelFromList(basemodellist, im)
    call mp%model_da()
    deallocate(mp)
  enddo
  !
  ! -- Deallocate for each exchange
  do ic = 1, baseexchangelist%Count()
    ep => GetBaseExchangeFromList(baseexchangelist, ic)
    call ep%exg_da()
    deallocate(ep)
  enddo
  !
  ! -- Deallocate for each solution
  do is=1,basesolutionlist%Count()
    sp => GetBaseSolutionFromList(basesolutionlist, is)
    call sp%sln_da()
    deallocate(sp)
  enddo
  !
  ! -- Deallocate solution group and simulation variables
  do isg = 1, solutiongrouplist%Count()
    sgp => GetSolutionGroupFromList(solutiongrouplist, isg)
    call sgp%sgp_da()
    deallocate(sgp)
  enddo
  call simulation_da()
  call lists_da()
  !
  ! -- Calculate memory usage, elapsed time and terminate
  call mem_usage(iout)
  call mem_da()
  call elapsed_time(iout, 1)
  call final_message()
  !
end program mf6

