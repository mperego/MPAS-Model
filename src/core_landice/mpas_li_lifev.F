!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_lifev
!
!> \MPAS land-ice velocity driver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id:$
!> \details
!>  This module contains the routines for interfacing with the 
!>  external velocity solver LifeV and optionally the PHG
!>  Full Stokes velocity solver that requires LifeV
!>
!
!-----------------------------------------------------------------------

module li_lifev


   use mpas_grid_types
   use mpas_configure
   use mpas_dmpar
   !use, intrinsic :: iso_c_binding

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: li_lifev_init, &
             li_lifev_block_init, &
             li_lifev_solve, &
             li_lifev_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************



contains



!***********************************************************************
!
!  routine li_lifev_init
!
!> \brief   Initializes velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the ice velocity solver in the LifeV
!>  external velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_lifev_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      integer :: err_tmp

      err = 0
      err_tmp = 0

      ! Check for configuration options that are incompatible with LifeV velocity solver
      if (config_num_halos /= 3) then
         write(0,*) "Error: LifeV velocity solvers requires that config_num_halos=3"
         err_tmp = 1
      endif
      err = ior(err,err_tmp)
      if (config_number_of_blocks /= 0) then
         write(0,*) "Error: LifeV velocity solvers requires that config_number_of_blocks=0"
         err_tmp = 1
      endif
      err = ior(err,err_tmp)

#ifdef USE_LIFEV
      ! These calls are needed for using any of the LifeV velocity solvers

      !call lifeV and set the grid of the velocity solver
      call velocity_solver_init_mpi(domain % dminfo % comm)

      if (config_dycore == 'phg-Stokes')
          call phg_init(domain, err)
      endif
      err = ior(err,err_tmp)


#else

      err = 1
      write(0,*) "Error: Compiling with LifeV required to run LifeV-L1L2, LifeV-FO and PHG-Stokes models"

#endif


      ! === error check
      if (err > 0) then
          write (0,*) "An error has occurred in li_lifev_init."
      endif

   !--------------------------------------------------------------------
   end subroutine li_lifev_init



!***********************************************************************
!
!  routine li_lifev_block_init
!
!> \brief   Initializes velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details 
!>  This routine initializes each block of the ice velocity solver in the 
!>  LifeV external velocity solver.
!
!-----------------------------------------------------------------------

   subroutine li_lifev_block_init(block, err)

       use mpas_timer

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (block_type), intent(in) :: &
         block          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: nCells, nEdges, nVertices, nCellsSolve, nEdgesSolve, nVerticesSolve, nVertLevels

      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnCell, verticesOnEdge

      integer, dimension(:), pointer :: indexToCellID

      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, areaTriangle

      ! halo exchange arrays
      integer, dimension(:), pointer :: sendCellsArray, &
                                     recvCellsArray, &
                                     sendVerticesArray, &
                                     recvVerticesArray, &
                                     sendEdgesArray, &
                                     recvEdgesArray

      err = 0

      ! Note: LifeV only supports one block per processor, but this has (hopefully)
      ! been written to work if that were to change.  (That's why all these LifeV init
      ! calls are in li_lifev_block_init instead of li_lifev_init.)

      !extract data from domain
      nCellsSolve = block % mesh % nCellsSolve
      nEdgesSolve = block % mesh % nEdgesSolve
      nVerticesSolve = block % mesh % nVerticesSolve

      nCells = block % mesh % nCells
      nEdges = block % mesh % nEdges
      nVertices = block % mesh % nVertices
      nVertLevels = block % mesh % nVertLevels

      cellsOnEdge => block % mesh % cellsOnEdge % array
      cellsOnVertex => block % mesh % cellsOnVertex % array
      verticesOnCell => block % mesh % verticesOnCell % array
      verticesOnEdge => block % mesh % verticesOnEdge % array
      indexToCellID => block % mesh % indexToCellID % array

      xCell => block % mesh % xCell % array
      yCell => block % mesh % yCell % array
      zCell => block % mesh % zCell % array

      areaTriangle => block % mesh % areaTriangle % array


      !build send and receive arrays using exchange_list
!      call array_from_exchange_list(block % parinfo % verticesToSend, sendVerticesArray)
!      call array_from_exchange_list(block % parinfo % verticesToRecv, recvVerticesArray)
!      call array_from_exchange_list(block % parinfo % cellsToSend, sendCellsArray)
!      call array_from_exchange_list(block % parinfo % cellsToRecv, recvCellsArray)
!      call array_from_exchange_list(block % parinfo % edgesToSend, sendEdgesArray)
!      call array_from_exchange_list(block % parinfo % edgesToRecv, recvEdgesArray)

#ifdef USE_LIFEV
      ! These calls are needed for using any of the LifeV velocity solvers

      !zCell is supposed to be zero for L1L2 and FO solvers
      !nVertLevels should be equal to nVertLevelsSolve (no split the domain in the vertical direction)
      call mpas_timer_start("velocity_solver_set_grid_data")
      call velocity_solver_set_grid_data(nCells, nEdges, nVertices, nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, &
	cellsOnEdge, cellsOnVertex, verticesOnCell, verticesOnEdge, indexToCellID, xCell, yCell, zCell, areaTriangle, &
	sendCellsArray, recvCellsArray, sendEdgesArray, recvEdgesArray, sendVerticesArray, recvVerticesArray)
      call mpas_timer_stop("velocity_solver_set_grid_data")
#else
      write(*,*) "LifeV needed to run L1L2, FO and Stokes models"
      write(0,*) "LifeV needed to run L1L2, FO and Stokes models"
      err = 1
      return
#endif

      !these can be deallocated because they have been copied on the c++ side
      deallocate(sendCellsArray, &
                 recvCellsArray, &
                 sendVerticesArray, &
                 recvVerticesArray, &
                 sendEdgesArray, &
                 recvEdgesArray)


   !--------------------------------------------------------------------
   end subroutine li_lifev_block_init



!***********************************************************************
!
!  routine li_lifev_solve
!
!> \brief   Interface to call LifeV velocity solvers
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details 
!>  This routine calls external LifeV velocity solvers and/or PHG velocity solvers.
!
!-----------------------------------------------------------------------

   subroutine li_lifev_solve(mesh, state, err)

       use mpas_timer

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (state_type), intent(inout) :: &
         state          !< Input: state information 

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      real (kind=RKIND), dimension(:), pointer :: &
         thickness, lowerSurface, upperSurface, beta, LayerThicknessFractions
      real (kind=RKIND), dimension(:,:), pointer :: &
         normalVelocity
      real (kind=RKIND), dimension(:,:,:), pointer :: &
         tracers
      integer, dimension(:), pointer :: vertexMask

      integer :: index_temperature
      real (kind=RKIND) :: anyVertexMaskChanged  ! THIS TYPE NEEDS TO BE CHANGE TO INTEGER ONCE REGISTRY IS UPDATED

      err = 0
      
      ! Mesh variables
      LayerThicknessFractions => mesh % LayerThicknessFractions % array
      beta => mesh % beta % array

      ! State variables
      normalVelocity => state % normalVelocity % array
      thickness => state % thickness % array
      lowerSurface => state % lowerSurface % array
      upperSurface => state % upperSurface % array
      tracers => state % tracers % array
      index_temperature = state % index_temperature
      vertexMask => state % vertexMask % array
      anyVertexMaskChanged = state % anyVertexMaskChanged % scalar


#ifdef USE_LIFEV
      ! LifeV calls to be made only when vertex mask changes
      if (anyVertexMaskChanged .eq. 1) then
          call mpas_timer_start("velocity_solver_compute_2d_grid")
          call velocity_solver_compute_2d_grid(vertexMask)
          call mpas_timer_stop("velocity_solver_compute_2d_grid")

          select case (config_dycore)
          case ('L1L2')
              call mpas_timer_start("velocity_solver_init_L1L2")
	      !call velocity_solver_extrude_3d_grid(LayerThicknessFractions, lowerSurface, thickness)
	      call velocity_solver_init_L1L2(LayerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_L1L2")
          case ('FO')
              call mpas_timer_start("velocity_solver_extrude_3d_grid")
              call velocity_solver_extrude_3d_grid(LayerThicknessFractions, lowerSurface, thickness)
              call mpas_timer_stop("velocity_solver_extrude_3d_grid")
              call mpas_timer_start("velocity_solver_init_FO")
              call velocity_solver_init_FO(LayerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_FO")
          case ('Stokes')
#ifdef USE_PHG
              call mpas_timer_start("velocity_solver_extrude_3d_grid")
              call velocity_solver_extrude_3d_grid(LayerThicknessFractions, lowerSurface, thickness)
              call mpas_timer_stop("velocity_solver_extrude_3d_grid")
              call mpas_timer_start("velocity_solver_init_stokes")
              call velocity_solver_init_stokes(LayerThicknessFractions)
              call mpas_timer_stop("velocity_solver_init_stokes")
#else
              write(0,*) "Error: PHG library needed to run stokes"
              err = 1
              return
#endif
          end select
      endif

      ! LifeV calls to be made every time step (solve velocity!)
      select case (config_dycore)
      case ('L1L2')
          call mpas_timer_start("velocity_solver_solve_L1L2")
          call velocity_solver_solve_L1L2(lowerSurface, thickness, beta, tracers(index_temperature,:,:), normalVelocity)
!         call velocity_solver_estimate_SS_SMB(normalVelocity, mesh % sfcMassBal % array)
          call mpas_timer_stop("velocity_solver_solve_L1L2")
          ! Optional calls to have LifeV output data files
          call mpas_timer_start("velocity_solver export")
          call velocity_solver_export_2d_data(lowerSurface, thickness, beta)
          call velocity_solver_export_L1L2_velocity();
          call mpas_timer_stop("velocity_solver export")
      case ('FO')
          call mpas_timer_start("velocity_solver_solve_FO")
          call velocity_solver_solve_FO(lowerSurface, thickness, beta, tracers(index_temperature,:,:), normalVelocity)
!	  call velocity_solver_estimate_SS_SMB(normalVelocity, mesh % sfcMassBal % array)
          call mpas_timer_stop("velocity_solver_solve_FO")
          call mpas_timer_start("velocity_solver export")
!	  call velocity_solver_init_L1L2(LayerThicknessFractions)          
!	  call velocity_solver_export_2d_data(lowerSurface, thickness, beta)
          call velocity_solver_export_FO_velocity()
          call mpas_timer_stop("velocity_solver export")
      case ('Stokes')
#ifdef USE_PHG
          call mpas_timer_start("velocity_solver_solve_stokes")
          call velocity_solver_solve_stokes(lowerSurface, thickness, beta, tracers(index_temperature,:,:), normalVelocity)
          call mpas_timer_stop("velocity_solver_solve_stokes")
#else
          write(0,*) "Error: PHG library needed to run stokes"
          err = 1
          return
#endif
      end select
#else
      write(0,*) "Error: LifeV needed to run L1L2, FO and Stokes models"
      err = 1
      return
#endif

   !--------------------------------------------------------------------
   end subroutine li_lifev_solve



!***********************************************************************
!
!  routine li_lifev_finalize
!
!> \brief   Initializes velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine finalizes the ice velocity solver in the LifeV external library.
!
!-----------------------------------------------------------------------

   subroutine li_lifev_finalize(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

#ifdef USE_LIFEV
      ! This call is needed for using any of the LifeV velocity solvers
      call velocity_solver_finalize()
#else
      write(*,*) "LifeV needed to run L1L2, FO and Stokes models"
      err = 1
      return
#endif

   !--------------------------------------------------------------------
   end subroutine li_lifev_finalize



!***********************************************************************
!  private subroutines
!***********************************************************************



!***********************************************************************
!
!  routine phg_init
!
!> \brief   Initializes phg velocity solver
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the ice velocity solver in the PHG
!>  external library.
!
!-----------------------------------------------------------------------

   subroutine phg_init(domain, err)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

#ifdef USE_PHG
      ! This call is needed for using any of the PHG velocity solvers
      call phg_init(domain % dminfo % comm)
#else
      write(0,*) "Error: Compiling with PHG required to run PHG-Stokes model"
      err = 1
      return
#endif

   !--------------------------------------------------------------------
   end subroutine phg_init



!***********************************************************************
!
!  routine array_from_exchange_list
!
!> \brief   Converts the MPAS Exchange Lists to flat arrays for external use
!> \author Matt Hoffman
!> \date   3 October 2013
!> \version SVN:$Id$
!> \details
!>  This routine converts the MPAS Exchange Lists to flat arrays for 
!>  use by external dycores.
!-----------------------------------------------------------------------

   subroutine array_from_exchange_list(exlist, array)

     implicit none

     type (exchange_list), pointer, intent(in) :: exlist
     type (exchange_list), pointer :: listPtr

     integer, dimension(:), pointer, intent(out) :: array
     integer :: dataSize, offset, i

!     dataSize = 1 !in first position we will store the size of the vector
!     listPtr => exlist
!     do while (associated(listPtr))
!       dataSize = dataSize + (listPtr % nlist) + 2
!       listPtr => listPtr % next
!     end do

!     allocate(array(dataSize))

!     array(1) = dataSize;
!     offset = 2;
!     listPtr => exlist
!     do while (associated(listPtr))
!       array(offset) = listPtr % procID
!       offset = offset + 1
!       array(offset) = listPtr % nlist
!       do i=1,listPtr % nlist
!         array(i+offset) = listPtr % list(i) -1
!       end do
!       offset = offset + listPtr % nlist + 1
!       listPtr => listPtr % next
!     end do

   end subroutine array_from_exchange_list

!***********************************************************************

end module li_lifev

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
