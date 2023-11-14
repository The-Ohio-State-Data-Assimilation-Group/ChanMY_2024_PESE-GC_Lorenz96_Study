! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! Module to implement Probit-space Ensemble Size Expansion (PESE; "peace").
! Written by Man-Yau ("Joseph") Chan

module pese_mod

   use types_mod,                only : r8, i4, i8
   
   use io_filenames_mod,         only : file_info_type, NO_IO
   
   use random_seq_mod,           only : random_seq_type, random_gaussian, init_random_seq

   use utilities_mod,            only : error_handler, E_ERR, E_MSG, E_DBG

   use ensemble_manager_mod,     only : ensemble_type

   use probit_transform_mod,     only : transform_to_probit, transform_from_probit

   use algorithm_info_mod,       only : probit_dist_info

   use location_mod,             only : location_type

   use distribution_params_mod,  only : distribution_params_type, NORMAL_DISTRIBUTION

   use assim_model_mod,          only : get_state_meta_data

   
   
   implicit none
   private
   
   public :: init_pese_info, generate_pese_resampling_weights, pese_type, &
             pese_modify_file_info_io_flags, apply_pese, pese_modify_all_file_info_io_flags
   

   ! Useful constants
   real(r8), parameter                    :: MEAN_ZERO = 0.
   real(r8), parameter                    :: SPREAD_UNITY = 1.

   ! Message string
   character(len=512)      :: msgstring
   character(len=*), parameter :: source = 'pese_mod.f90'

   
   ! Structured type holding information needed for PESE
   type pese_type
      
      ! Number of ensemble members produced by time-integrating forecast model
      integer                             :: ens_size_dynamic
      
      ! Number of PESE-created ensemble members
      integer                             :: ens_size_virtual
   
      ! Expanded ensemble size
      integer                             :: ens_size_expanded
   
      ! Resampling weights used to create PESE members
      real,                allocatable    :: resampling_weights(:,:)
   
      ! Matrices used in constructing resampling weights
      real,                allocatable    :: WWT      (:,:)
      real,                allocatable    :: L_W      (:,:)
      real,                allocatable    :: inv_L_W  (:,:)
      real,                allocatable    :: L_E      (:,:)
      real,                allocatable    :: avg_of_noise_columns             (:)
   
      ! k factor reference in Chan Anderson Chen 2020
      real                                :: k_factor

      ! Arrays to hold univariate ensmebles 
      real(r8),            allocatable    :: ens_dyn_1d(:)
      real(r8),            allocatable    :: ens_vir_1d(:)

      ! Variable controlling random number generator
      type(random_seq_type)               :: rng_info

      ! Flag to indicate if pese is used
      logical                             :: use_pese = .false.

   
   end type pese_type
   

   ! Allocate pointer for pese_type
   type( pese_type )       :: pese_info
   
   
   
   
   contains







   ! Subroutine to apply PESE resampling
   subroutine apply_pese( ens_handle, print_pese_debug )

      type(ensemble_type),          intent(inout)  :: ens_handle
      logical,                      intent(in)     :: print_pese_debug
      integer                                      :: my_state_kind
      type(location_type)                          :: my_state_loc
      integer                                      :: v_ind, imem
      logical                                      :: bounded_below, bounded_above
      real(r8)                                     :: lower_bound,   upper_bound
      integer                                      :: dist_type
      type(distribution_params_type), allocatable  :: dist_params
      type(distribution_params_type), allocatable  :: dist_params_copies(:)
      real(r8)                                     :: probit_avg1, probit_std1
      real(r8)                                     :: probit_avg2, probit_std2
      integer                                      :: num_ens_loop, iloop
      real(r8)                                     :: ratio_virtual_to_dynamic
      integer                                      :: ie_st, ie_ed

      real(r8)                                     :: original_cov, pese_cov


      ! TO DO LIST:
      ! 1) Incorporate "groups" into PESE




      ! Compute resampling coefficients
      pese_info%resampling_weights = 0.
      call generate_pese_resampling_weights()

      
      ! Loop over all variables in ensemble handle
      PESE_LOOP_OVER_VARIABLES: do v_ind = 1, ens_handle%my_num_vars

         allocate( dist_params )

         ! Obtain meta data for current state element
         call get_state_meta_data( ens_handle%my_vars( v_ind ), my_state_loc, my_state_kind ) 

         ! print out pre-probit avg and spread 
         if (print_pese_debug .and. v_ind == 1) then
            probit_avg1 = sum( ens_handle%copies(1:pese_info%ens_size_dynamic, v_ind ) ) / pese_info%ens_size_dynamic
            probit_std1 = sum( ( ens_handle%copies(1:pese_info%ens_size_dynamic, v_ind ) - probit_avg1 )**2 ) / ( pese_info%ens_size_dynamic - 1. )
            probit_std1 = sqrt( probit_std1 )
            write(*,'(x,a80,x,f,x,f)') 'Pre-PESE avg and standard dev', probit_avg1, probit_std1
         end if 
         

         ! Load user-selected marginal
         call probit_dist_info( my_state_kind, .true., .false., dist_type, &
               bounded_below, bounded_above, lower_bound, upper_bound )   

         ! Transform to probit space.
         call transform_to_probit(                                                                 &
            pese_info%ens_size_dynamic, ens_handle%copies(1:pese_info%ens_size_dynamic, v_ind ),   &
            dist_type, dist_params, pese_info%ens_dyn_1d(:), .false.,                              &
            bounded_below, bounded_above, lower_bound, upper_bound                                 &
         )
         
         ! TEMPORARY ADDITION START
         ! ------------------------
         ! if ( v_ind == 1 ) then
         !    write(*,*) 'Probit-transformed'
         !    write(*,*) pese_info%ens_dyn_1d(:)
         !    write(*,*) 'original'
         !    write(*,*) ens_handle%copies(1:pese_info%ens_size_dynamic, v_ind )
         ! end if
         !pese_info%ens_dyn_1d = ens_handle%copies(1:pese_info%ens_size_dynamic, v_ind )
         ! TEMPORARY ADDITION END
         ! ----------------------

         ! Check probit-transformed ensemble mean and spread.
         probit_avg1 = sum( pese_info%ens_dyn_1d(:) ) / pese_info%ens_size_dynamic
         probit_std1 = sum( ( pese_info%ens_dyn_1d(:) - probit_avg1 )**2 ) / ( pese_info%ens_size_dynamic - 1. )
         probit_std1 = sqrt( probit_std1 )
         if (print_pese_debug .and. v_ind == 1) then
            write(*,'(x,a80,x,f6.2,x,f6.2)') 'Probit-transformed ensemble avg and standard deviation', probit_avg1, probit_std1
         end if

         ! Ensure probit-transformed ensemble mean is zero and spread is unity
         REGULARIZE_PROBIT_STATS: if ( dist_type .ne. NORMAL_DISTRIBUTION ) then
            pese_info%ens_dyn_1d(:) = pese_info%ens_dyn_1d(:) - probit_avg1
            pese_info%ens_dyn_1d(:) = pese_info%ens_dyn_1d(:)/probit_std1
         end if REGULARIZE_PROBIT_STATS

         ! Apply PESE assuming copula is Gaussian
         call apply_pese_for_gaussian_copula()

         ! Check post-PESE avg and std dev
         if (print_pese_debug .and. v_ind == 1) then
            probit_avg2 = ( sum( pese_info%ens_dyn_1d(:) ) + sum( pese_info%ens_vir_1d(:) ) ) / pese_info%ens_size_expanded
            probit_std2 = sum( ( pese_info%ens_dyn_1d(:) - probit_avg2 )**2 ) + sum( ( pese_info%ens_vir_1d(:) - probit_avg2 )**2 )
            probit_std2 = sqrt( probit_std2 / ( pese_info%ens_size_expanded-1. ) )
            write(*,'(x,a80,x,f6.2,x,f6.2)') 'PESE-expanded probit-transformed ensemble avg and standard deviation', probit_avg2, probit_std2
         end if

         ! Because inverse probit transform subroutine can only accept ens_size_dynamic probit values, to process ens_size_virtual probit 
         ! values, need to run transform_from_probit multiple times
         ratio_virtual_to_dynamic = pese_info%ens_size_virtual/ pese_info%ens_size_dynamic
         num_ens_loop = ceiling( ratio_virtual_to_dynamic )

         ! Make many copies of dist_params
         allocate( dist_params_copies(num_ens_loop+2) )
         dist_params_copies(:) = dist_params


         ! Invert PPI transform for all but last few virtual members
         INV_PPI_FOR_MOST_VIRTUAL_MEMBERS: do iloop = 1, num_ens_loop-1

            ! Determine start and end indices of virtual members to process
            ie_st = (iloop-1) * pese_info%ens_size_dynamic
            ie_ed =     iloop * pese_info%ens_size_dynamic    

            ! Invert PPI for selected virtual members
            call transform_from_probit(                                                                                 &
               pese_info%ens_size_dynamic, pese_info%ens_vir_1d( (ie_st+1):ie_ed  ), dist_params_copies(iloop),             &
               ens_handle%copies( (pese_info%ens_size_dynamic+ie_st+1):(pese_info%ens_size_dynamic+ie_ed+1), v_ind )    &
            )
            

         end do INV_PPI_FOR_MOST_VIRTUAL_MEMBERS


         ! Invert PPI for the last ens_size_dynamic members
         ie_st = pese_info%ens_size_virtual - pese_info%ens_size_dynamic
         ie_ed = pese_info%ens_size_virtual
         call transform_from_probit(                                                                                 &
               pese_info%ens_size_dynamic, pese_info%ens_vir_1d( (ie_st+1):ie_ed  ), dist_params_copies(num_ens_loop),                           &
               ens_handle%copies( (pese_info%ens_size_dynamic+ie_st+1):(pese_info%ens_size_dynamic+ie_ed+1), v_ind )    &
         )


            

         ! TEMPORARY ADDITION START
         ! ------------------------
         ! if (print_pese_debug .and. v_ind == 1) then
         !    write(*,*)'inv-PPI vs no transform'
         !       write(*,*) ens_handle%copies( pese_info%ens_size_dynamic+imem, v_ind ), pese_info%ens_vir_1d(imem)
         !    end do
         ! endif
         !ens_handle%copies( (pese_info%ens_size_dynamic+1):pese_info%ens_size_expanded, v_ind ) = pese_info%ens_vir_1d(:)
         ! TEMPORARY ADDITION END
         ! ----------------------

         
         ! print out post-probit avg and spread 
         if (print_pese_debug .and. v_ind == 1) then
            probit_avg1 = sum( ens_handle%copies(1:pese_info%ens_size_expanded, v_ind ) ) / pese_info%ens_size_expanded
            probit_std1 = sum( ( ens_handle%copies(1:pese_info%ens_size_expanded, v_ind ) - probit_avg1 )**2 ) 
            probit_std1 = probit_std1 / ( pese_info%ens_size_expanded - 1. )
            probit_std1 = sqrt( probit_std1 )
            write(*,'(x,a80,x,f,x,f)') 'Post-PESE avg and standard dev', probit_avg1, probit_std1
         end if

         
         ! Deallocate dist_params in preparation for the next iteration of the PESE_LOOP_OVER_VARIABLES
         deallocate( dist_params_copies )
         deallocate( dist_params )        



      end do PESE_LOOP_OVER_VARIABLES  




      ! Compute ens covariances before and after PESE
      if (print_pese_debug ) then

         ! Original cov
         probit_avg1 = sum( ens_handle%copies(1:pese_info%ens_size_dynamic, 1 ) ) / pese_info%ens_size_dynamic
         probit_avg2 = sum( ens_handle%copies(1:pese_info%ens_size_dynamic, 2 ) ) / pese_info%ens_size_dynamic
         original_cov =  sum(                                                             &
            (ens_handle%copies(1:pese_info%ens_size_dynamic, 1 ) - probit_avg1)           &
            * (ens_handle%copies(1:pese_info%ens_size_dynamic, 2 ) - probit_avg2)         &
         ) / ( pese_info%ens_size_dynamic -1. )

         ! New cov
         probit_avg1 = sum( ens_handle%copies(1:pese_info%ens_size_expanded, 1 ) ) / pese_info%ens_size_expanded
         probit_avg2 = sum( ens_handle%copies(1:pese_info%ens_size_expanded, 2 ) ) / pese_info%ens_size_expanded
         pese_cov =  sum(                                                             &
            (ens_handle%copies(1:pese_info%ens_size_expanded, 1 ) - probit_avg1)           &
            * (ens_handle%copies(1:pese_info%ens_size_expanded, 2 ) - probit_avg2)         &
         ) / ( pese_info%ens_size_expanded -1. )

         write(*,'(x,a80,x,f)') 'Pre-PESE covariance', original_cov
         write(*,'(x,a80,x,f)') 'Aft-PESE covariance', pese_cov

      end if


   end subroutine apply_pese

















   ! Subroutine to apply PESE resampling strategy for Gaussian marginals
   ! NOTE: THIS SUBROUTINE ONLY RESAMPLES A UNIVARIATE ENSEMBLE INPUTTED INTO PESE_INFO
   subroutine apply_pese_for_gaussian_copula( )

      ! Indexing variables
      integer                                :: i,j,k

      ! Variables to store mean and spread
      real(r8)                               :: var_dyn, var_exp    ! Variances of dynamic and expanded ensembles
      real(r8)                               :: avg_dyn, avg_vir    ! Avgs of dynamic and virtual ensembles


      ! Convert dynamic ensemble to dynamic perturbations
      avg_dyn = sum( pese_info%ens_dyn_1d ) / pese_info%ens_size_dynamic
      pese_info%ens_dyn_1d = pese_info%ens_dyn_1d - avg_dyn

      ! Apply resampling coefficients to generate virtual perturbations
      GAUSSIAN_RESAMPLING: do i = 1, pese_info%ens_size_virtual
         pese_info%ens_vir_1d(i) = sum( pese_info%ens_dyn_1d         &
                                        * pese_info%resampling_weights(:,i) )
      end do GAUSSIAN_RESAMPLING

      ! Check pre-pese and post-pese avgs
      avg_vir = sum( pese_info%ens_vir_1d ) / pese_info%ens_size_virtual
      if ( abs(avg_vir) > abs(avg_dyn)/1e5 .and. abs(avg_vir) > 1e-5 ) then
         write(*,*) '       Post-PESE perturbation avg value: ', avg_vir
         write(msgstring,*) 'ERROR: Post-PESE perturbations have non-zero mean!'
         call error_handler(E_ERR,'apply_pese_for_gaussian_marginal', msgstring, source)
      end if

      ! Check pre-pese and post-pese variances
      var_dyn = sum( pese_info%ens_dyn_1d**2 )
      var_exp = var_dyn + sum( pese_info%ens_vir_1d**2 )
      var_dyn = var_dyn / (pese_info%ens_size_dynamic -1.)
      var_exp = var_exp / (pese_info%ens_size_expanded-1.)
      if ( abs(var_exp - var_dyn) > var_dyn/1e5 ) then
         write(*,*) '       Post-PESE variance: ', var_exp
         write(*,*) '        Pre-PESE variance: ', var_dyn          
         write(msgstring,*) 'ERROR: Post-PESE ensemble variance does not match pre-PESE ensemble variance!'
         call error_handler(E_ERR,'apply_pese_for_gaussian_marginal', msgstring, source)
      end if

      ! Convert perturbations to ensemble
      pese_info%ens_dyn_1d = pese_info%ens_dyn_1d + avg_dyn
      pese_info%ens_vir_1d = pese_info%ens_vir_1d + avg_dyn

   end subroutine apply_pese_for_gaussian_copula
   
















   
   
   ! Subroutine to initialize variables needed for PESE
   subroutine init_pese_info( expanded_ens_size, dynamic_ens_size )
   
      integer,                intent(in)        :: expanded_ens_size, dynamic_ens_size
      integer                                   :: virtual_ens_size
   
      ! Compute number of members created by PESE
      virtual_ens_size = expanded_ens_size - dynamic_ens_size
   
      ! Allocate useful arrays for resampling weight calculations
      allocate( pese_info%resampling_weights      ( dynamic_ens_size, virtual_ens_size ) )
      allocate( pese_info%WWT                     ( dynamic_ens_size, dynamic_ens_size ) )  ! Aka, matrix C_w
      allocate( pese_info%L_W                     ( dynamic_ens_size, dynamic_ens_size ) )  !
      allocate( pese_info%inv_L_W                 ( dynamic_ens_size, dynamic_ens_size ) )  ! Aka, matrix (C_w)^-0.5
      allocate( pese_info%L_E                     ( dynamic_ens_size, dynamic_ens_size ) )  ! Aka, matrix (C_E)^+0.5
      allocate( pese_info%avg_of_noise_columns    ( dynamic_ens_size )                   )
   
      ! Fill with zeros
      pese_info%resampling_weights               = 0.
      pese_info%WWT    = 0.
      pese_info%L_W = 0.
      pese_info%inv_L_W = 0.
      pese_info%L_E             = 0.
   
      ! Fill in ensemble sizes
      pese_info%ens_size_dynamic    =  dynamic_ens_size
      pese_info%ens_size_virtual    =  virtual_ens_size
      pese_info%ens_size_expanded   =  expanded_ens_size
   
      ! Fill in recaling factor k
      pese_info%k_factor            = sqrt( (expanded_ens_size-1.) / (dynamic_ens_size-1.) )

      ! Allocate arrays to hold univariate ensembles
      allocate( pese_info%ens_dyn_1d   (dynamic_ens_size) )
      allocate( pese_info%ens_vir_1d   (virtual_ens_size) )

      ! Initialize random number generator 
      call init_random_seq( pese_info%rng_info, 9999 )

      ! Indicate that PESE is being used
      pese_info%use_pese=.true.
   
   end subroutine init_pese_info
   

   
   ! Subroutine to modify IO copy flags for 6 kinds of file_info
   subroutine pese_modify_all_file_info_io_flags(dynamic_ens_size, ens_size, ENS_MEM_START, file_info_input, file_info_forecast, file_info_preassim, file_info_postassim, file_info_analysis, file_info_output)

      integer,                            intent(in)     :: dynamic_ens_size
      integer,                            intent(in)     :: ens_size
      integer,                            intent(in)     :: ENS_MEM_START
      type( file_info_type ),             intent(inout)  :: file_info_input
      type( file_info_type ),             intent(inout)  :: file_info_forecast
      type( file_info_type ),             intent(inout)  :: file_info_preassim
      type( file_info_type ),             intent(inout)  :: file_info_postassim
      type( file_info_type ),             intent(inout)  :: file_info_analysis
      type( file_info_type ),             intent(inout)  :: file_info_output

      call pese_modify_file_info_io_flags( file_info_input, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )
   
      call pese_modify_file_info_io_flags( file_info_forecast, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )
   
      call pese_modify_file_info_io_flags( file_info_preassim, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )
   
      call pese_modify_file_info_io_flags( file_info_postassim, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )
   
      call pese_modify_file_info_io_flags( file_info_analysis, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )
   
      call pese_modify_file_info_io_flags( file_info_output, ENS_MEM_START, &
                                          dynamic_ens_size, ens_size )

   end subroutine pese_modify_all_file_info_io_flags

   









   
   
   
   
   ! Subroutine to overwrite IO copy flags in file_info
   ! --------------------------------------------------
   ! Goal is to adjust read_state() and write_state() to ignore members that are generated by PESE (virtual members)
   subroutine pese_modify_file_info_io_flags( file_info, ENS_MEM_START, ens_size_dynamic, ens_size_expanded )

      type(file_info_type),         intent(inout)  :: file_info
      integer,                      intent(in)     :: ENS_MEM_START
      integer,                      intent(in)     :: ens_size_dynamic
      integer,                      intent(in)     :: ens_size_expanded
      integer                                      :: c0, c1
      
      ! Overwrite IO flags for PESE-generated members
      c0 = ENS_MEM_START + ens_size_dynamic
      c1 = ENS_MEM_START + ens_size_expanded -1
      file_info%stage_metadata%io_flag( c0:c1 ) = NO_IO

   end subroutine pese_modify_file_info_io_flags
   
   
   
   
   

   
   
   ! Subroutine to generate probit-space resampling weights for virtual member creation
   ! ----------------------------------------------------------------------------------
   ! These weights are s.t. the expanded ensemble (dynamical+virtual members) has the same 
   ! Gaussian copula as the original ensemble.
   ! Algorithm is described in Appendix B of Chan, Anderson, Chen (2020; MWR)
   ! Side-note: virtual members are also known as "PESE-created" members
   subroutine generate_pese_resampling_weights(  )
   
      real                                         :: weight_offset        ! Aka, bar(E) in Chan, Anderson, Chen (2020; MWR)
      integer                                      :: i, j                 ! Dummy indices
      integer                                      :: error_flag
      real,                         allocatable    :: buffer2d(:,:)
   
      ! Draw Gaussian white noise
      do i = 1, pese_info%ens_size_dynamic
         do j = 1, pese_info%ens_size_virtual
            pese_info%resampling_weights( i, j )         &
               = random_gaussian( pese_info%rng_info, MEAN_ZERO, SPREAD_UNITY )
         end do
      end do
   
   
      ! Remove average of columns
      pese_info%avg_of_noise_columns &
         = sum( pese_info%resampling_weights, 2) / pese_info%ens_size_virtual
      do i = 1, pese_info%ens_size_dynamic
         pese_info%resampling_weights( i, : ) &
            = pese_info%resampling_weights( i, : ) - pese_info%avg_of_noise_columns(i)
      end do
      
   
      ! Compute white noise outer product
      do i = 1, pese_info%ens_size_dynamic
         do j = 1, pese_info%ens_size_dynamic
            pese_info%WWT(i,j) &
               = sum( pese_info%resampling_weights( i, : ) &
                        * pese_info%resampling_weights( j, : )  )
         end do
      end do
  
   
      ! Generate noise decorrelation operator 
      pese_info%L_W(:,:)  = pese_info%WWT(:,:)
   
   
      ! Call Cholesky decomposition and then zero unneeded elements
      call spotrf( 'L', pese_info%ens_size_dynamic,                                    &
                  pese_info%L_W,                                                       &
                  pese_info%ens_size_dynamic, error_flag                               &
                  )
      do i = 1, pese_info%ens_size_dynamic
         pese_info%L_W( i, (i+1):pese_info%ens_size_dynamic ) = 0.
      end do
   
   
      ! Generate decorrelation operator by inverting Chol( white noise covariance )
      pese_info%inv_L_W(:,:) = pese_info%L_W(:,:)
      call strtri( 'L', 'N', pese_info%ens_size_dynamic,                &
                  pese_info%inv_L_W,                                    & 
                  pese_info%ens_size_dynamic, error_flag )



      do i = 1, pese_info%ens_size_dynamic 
         do j=1, pese_info%ens_size_dynamic 
            pese_info%WWT(i,j) = sum( pese_info%inv_L_W(i,:) * pese_info%L_W(:,j))
         end do
      end do
   
   
      ! Decorrelate white noise
      allocate( buffer2d      ( pese_info%ens_size_dynamic, pese_info%ens_size_virtual )     )
      buffer2d=0
      do i = 1, pese_info%ens_size_dynamic
         do j = 1, pese_info%ens_size_virtual
            buffer2d(i,j) = sum(                                  &
               pese_info%inv_L_W(i,:)                             &
               * pese_info%resampling_weights(:,j)                &
            )
         end do
      end do
      pese_info%resampling_weights(:,:) = buffer2d(:,:)
      
   
      ! CHECK: white noise covariance matrix
      do i = 1, pese_info%ens_size_dynamic
         do j = 1, pese_info%ens_size_dynamic
            pese_info%WWT(i,j) &
               = sum( pese_info%resampling_weights( i, : )        &
                        * pese_info%resampling_weights( j, : )  )
         end do
      end do



      ! Generate covariance that resampling weights must satisfy
      pese_info%L_E(:,:)                            &
         =  -1* (pese_info%k_factor-1.)**2 / (pese_info%ens_size_virtual)
   
      do i = 1, pese_info%ens_size_dynamic
         pese_info%L_E(i,i)                                                            &
            = ( pese_info%L_E(i,i)                                                     &
               + pese_info%ens_size_virtual / (pese_info%ens_size_dynamic-1.) )
      end do
   
   
      ! Generate correlating operator for resampling weights
      call spotrf( 'L', pese_info%ens_size_dynamic,                                    &
                  pese_info%L_E,                                                       &
                  pese_info%ens_size_dynamic, error_flag                               &
      )
      do i = 1, pese_info%ens_size_dynamic
         pese_info%L_E( i, (i+1):pese_info%ens_size_dynamic ) = 0.
      end do
   
   
      ! Apply correlating operator on resampling weights
      buffer2d(:,:) = 0.
      do i = 1, pese_info%ens_size_dynamic
         do j = 1, pese_info%ens_size_virtual
            buffer2d(i,j) = sum(                                                 &
               pese_info%L_E(i,:)                                                &
               * pese_info%resampling_weights(:,j)                               &
            )
         end do
      end do
      pese_info%resampling_weights(:,:) = buffer2d(:,:)
   
   
      ! Offset resampling weights
      pese_info%resampling_weights(:,:)                                          &
         = pese_info%resampling_weights(:,:)                                     &
            + (pese_info%k_factor-1.)/(pese_info%ens_size_virtual)
   
     
   
   end subroutine generate_pese_resampling_weights
   





   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   !========================================================================
   ! end module pese_mod
   !========================================================================
   
   end module pese_mod
   
   