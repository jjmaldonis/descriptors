
! Reverse Monte Carlo structural simulation
!
! Performs RMC refinement against reduced density
! function data from electron, x-ray, or neutron
! diffraction data and against fluctuation electron
! microscopy V(k) data.
!
! Voyles research group: Jason Maldonis, Jinwoo Hwang, Feng, Yi, and
! Paul Voyles, begun 12/16/08
! 

program rmc

#ifdef USE_LMP
    use LAMMPS
#endif
    use omp_lib
    use rmc_global
    use readinputs
    use model_mod
    use rmc_functions
    use eam_mod
    use vor_mod
    implicit none
    include 'mpif.h'
    ! LAMMPS objects
#ifdef USE_LMP
    type (C_ptr) :: lmp
    real (C_double), pointer :: te1 => NULL()
    real (C_double), pointer :: te2 => NULL()
    character (len=512) :: lmp_cmd_str
#endif
    !integer :: prev_moved_atom = 0, n
    ! RMC / Femsim objects
    type(model) :: m
    character (len=256) :: model_filename
    character (len=256) :: paramfile
    character (len=256) :: eam_filename
    character (len=256) :: outbase
    character (len=256) :: jobID, c, step_str
    character (len=512) :: comment
    character (len=256) :: time_elapsed, vki_fn, vku_fn, vkf_fn, output_model_fn, energy_fn, final_model_fn, chi_squared_file, acceptance_rate_fn, beta_fn
    real :: temperature
    real :: max_move
    real :: Q, res, alpha
    real :: cutoff
    real, pointer, dimension(:) :: vk, vk_exp, k, vk_exp_err, v_background, vk_as !vk_as is for autoslice/multislice (whatever it's called)
    real, pointer, dimension(:,:) :: cutoff_r
    real, pointer, dimension(:,:) :: scatfact_e
    real :: xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new
    real :: scale_fac, scale_fac_initial, beta
    double precision :: chi2_old, chi2_new, del_chi, chi2_gr, chi2_vk, chi2_no_energy, chi2_initial
    real :: R
    integer :: i, j, jj, jjj, step_start
    integer :: w
    integer :: nk
    integer :: ntheta, nphi, npsi
    integer :: istat, status2, length
    integer :: iseed2
    real :: randnum
#ifndef USE_LMP
    double precision :: te1, te2 ! For eam, not lammps version of energy
#endif
    logical :: square_pixel, accepted, use_rmc, use_multislice
    integer :: ipvd, nthr
    doubleprecision :: t0, t1, t2 !timers
    real :: x! This is the parameter we will use to fit vsim to vas.
    integer, dimension(100) :: acceptance_array
    real :: avg_acceptance = 1.0
    integer , dimension(maxneighs*2) :: update_these

    !------------------- Program setup. -----------------!

#ifdef TIMING
    write(*,*) "Using timing version of code!"
#endif
#ifdef USE_LMP
    write(*,*) "Using LAMMPS!"
#endif

    call mpi_init_thread(MPI_THREAD_MULTIPLE, ipvd, mpierr) !http://www.open-mpi.org/doc/v1.5/man3/MPI_Init_thread.3.php
    call mpi_comm_rank(mpi_comm_world, myid, mpierr)
    call mpi_comm_size(mpi_comm_world, numprocs, mpierr)

    !call omp_set_num_threads(2)
    nthr = omp_get_max_threads()
    if(myid.eq.0)then
        write(*,*)
        write(*,*) "Using", numprocs, "processors."
        write(*,*) "OMP found a max number of threads of", nthr
        write(*,*)
    endif

#ifdef USE_LMP
    if(myid.eq.0) then
        call lammps_open ('lmp -log log.simple -screen none', mpi_comm_world, lmp)
        call lammps_file (lmp, 'lmp_energy.in')
    endif
#endif

if(myid .eq. 0) then
    call get_command_argument(1, c, length, istat)
    if (istat == 0) then
        jobID = "_"//trim(c)
    else
        error stop "istat for jobid get_command_arg was nonzero"
        !jobID = '_temp'
    end if
    call get_command_argument(2, c, length, istat)
    if (istat == 0) then
        paramfile = trim(c)
    else
        error stop "istat for paramfile get_command_arg was nonzero"
    end if
    paramfile = trim(paramfile)
    call get_command_argument(3, c, length, istat)
    if (istat == 0) then
        model_filename = trim(c)
    else
        error stop "istat for modelfile get_command_arg was nonzero"
    end if
    model_filename = trim(model_filename)
    write(*,*) "Modelfile:", model_filename
    write(*,*) "Paramfile:", paramfile
    write(*,*) "JobID:", jobID

    ! Set output filenames.
    outbase = ""
    write(time_elapsed, "(A12)") "time_elapsed"
    time_elapsed = trim(trim(time_elapsed)//jobID)//".txt"
    write(output_model_fn, "(A12)") "model_update"
    output_model_fn = trim(trim(output_model_fn)//jobID)//".txt"
    write(final_model_fn, "(A11)") "model_final"
    final_model_fn = trim(trim(final_model_fn)//jobID)//".txt"
    !write(energy_fn, "(A6)") "energy"
    !energy_fn = trim(trim(energy_fn)//jobID)//".txt"
    write(chi_squared_file, "(A11)") "chi_squared"
    chi_squared_file = trim(trim(chi_squared_file)//jobID)//".txt"
    write(acceptance_rate_fn, "(A15)") "acceptance_rate"
    acceptance_rate_fn = trim(trim(acceptance_rate_fn)//jobID)//".txt"
endif

    !------------------- Read inputs and initialize. -----------------!

    ! Set input filenames.
    !param_filename = 'param_file.in'
    call mpi_bcast(paramfile, 256, MPI_CHARACTER, 0, mpi_comm_world, mpierr)
    call mpi_bcast(model_filename, 256, MPI_CHARACTER, 0, mpi_comm_world, mpierr)
    
    ! Start timer.
    t0 = omp_get_wtime()

    ! Read input model
    !call read_model(model_filename, comment, m, istat)
    call read_model(model_filename, m, istat)
    call check_model(m, istat)
    call recenter_model(0.0, 0.0, 0.0, m)

    ! Read input parameters
    allocate(cutoff_r(m%nelements,m%nelements),stat=istat)
    cutoff_r = 2.272
    !call read_inputs(param_filename,model_filename, eam_filename, step_start, temperature, max_move, cutoff_r, alpha, vk_exp, k, vk_exp_err, v_background, ntheta, nphi, npsi, scale_fac_initial, Q, status2)
    temperature = 86.5 ! temperature*(sqrt(0.7)**(step_start/200000))
    max_move = 1.5 !max_move*(sqrt(0.94)**(step_start/200000))
    alpha = 0.00025

    if(myid .eq. 0) then
    write(*,*) "Model filename: ", trim(model_filename)
    write(*,*)
    endif

    use_rmc = .true.
    beta=1./((8.6171e-05)*temperature)
    iseed2 = 104756
    if(myid.eq.0) write(*,*) "random number generator seed =", iseed2

#ifdef USE_LMP
    call lammps_command (lmp, 'run 0')
    call lammps_extract_compute (te1, lmp, 'pot', 0, 0)
#else
    eam_filename = 'ZrCuAl2011.eam.alloy'
    call read_eam(m,eam_filename)
    call eam_initial(m,te1)
#endif
    te1 = te1/m%natoms
    if(myid .eq. 0) write(*,*) "Energy = ", te1


    ! Calculate VP for every atom
    cutoff = 3.6
    call vor_init(m, cutoff, paramfile)
    do i=1, m%natoms
        call vp_atom(m, i, cutoff)
    enddo
    call group_indexes
    call check_neighs(m,cutoff)

    t1 = omp_get_wtime()

    !------------------- Start RMC. -----------------!

    call mpi_barrier(mpi_comm_world, mpierr)

    if(use_rmc) then ! End here if we only want femsim. Set the variable above.

        ! Calculate initial chi2
        chi2_no_energy = chi_square(alpha,numtypes, numtypes_sim)

        chi2_initial = chi2_no_energy
        chi2_old = chi2_no_energy + te1
#ifndef USE_LMP
        e2 = e1 ! eam
#endif

        i=0!step_start
        if(myid.eq.0)then
            write(*,*)
            write(*,*) "Initialization complete. Starting Monte Carlo."
            write(*,*) "Initial Conditions:"
            write(*,*) "   Step =       ", i
            write(*,*) "   Energy =     ", te1
            write(*,*) "   LSqF V(k) =  ", chi2_no_energy
            write(*,*) "   Temperature =", temperature
            write(*,*) "   Max Move=", max_move
            write(*,*)
            ! Reset time_elapsed, energy_function, chi_squared_file
#ifdef TIMING
            open(35,file=trim(time_elapsed),form='formatted',status='unknown')
                t1 = omp_get_wtime()
                write(35,*) numprocs, "processors are being used."
                write(35,*) "Step, Time elapsed, Avg time per step, This step's time"
            close(35)
#endif
            open(36,file=trim(chi_squared_file),form='formatted',status='unknown')
                write(36,*) "step, chi2, energy"
                write(36,*) i, chi2_no_energy, te1
            close(36)
            open(37,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                write(37,*) "step, acceptance rate averaged over last 1000 steps"
            close(37)
        endif


        t0 = omp_get_wtime()
        ! RMC loop begins. The loop never stops.
        !do while (i .lt. step_start+400000)
        do while (i .ge. 0)
        !do while (i .lt. 2)
#ifdef TIMING
            t2 = omp_get_wtime()
#endif

            if(myid .eq. 0) write(*,*) "Starting step", i

#ifdef TIMING
            if( i > 100) then
                if(myid .eq. 0) write(*,*) "STOPPING MC AFTER 100 STEPS"
                call mpi_finalize(mpierr)
                stop ! Stop after 100 steps for timing runs.
            endif
#endif

            !call check_neighs(m,cutoff)
            call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            ! check_curoffs returns false if the new atom placement is too close to
            ! another atom. Returns true if the move is okay. (hard shere cutoff)
            do while( .not. check_cutoffs(m,cutoff_r,w) )
                ! Check_cutoffs returned false so reset positions and try again.
                m%xx%ind(w) = xx_cur
                m%yy%ind(w) = yy_cur
                m%zz%ind(w) = zz_cur
                call random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, max_move)
            end do
            ! Update hutches, data for chi2, and chi2/del_chi
            call hutch_move_atom(m,w,xx_new, yy_new, zz_new)
    
#ifdef USE_LMP
            write(lmp_cmd_str, "(A9, I4, A3, F, A3, F, A3, F)") "set atom ", w, " x ", xx_new, " y ", yy_new, " z ", zz_new
            call lammps_command(lmp, trim(lmp_cmd_str))
            call lammps_command (lmp, 'run 0')
            call lammps_extract_compute (te2, lmp, 'pot', 0, 0)
#else
            call eam_mc(m, w, xx_cur, yy_cur, zz_cur, xx_new, yy_new, zz_new, te2)
#endif
            te2 = te2/m%natoms
            !if(myid .eq. 0) write(*,*) "Energy = ", te2

            ! Calculate a randnum for accept/reject
            randnum = ran2(iseed2)
            ! Decide whether to reject just based on the energy
            accepted = .true.
            if(accepted) then
            !write(*,*) "Old VP index of moved atom:"
            !call print_index(w)
            update_these = -1
            update_these(1) = w
            do jj=1, nneighs(w)
                do j=1, maxneighs*2
                    if(update_these(j) == neighs(w,jj)) exit
                    if(update_these(j) == -1) then
                        update_these(j) = neighs(w,jj)
                        exit
                    endif
                enddo
            enddo
            do jj=1, nneighs(w)
                do jjj=1, nneighs(neighs(w,jj))
                    do j=1, maxneighs*2
                        if(update_these(j) == neighs(neighs(w,jj),jjj) ) exit
                        if(update_these(j) == -1) then
                            update_these(j) = neighs(neighs(w,jj),jjj)
                            exit
                        endif
                    enddo
                enddo
            enddo
            do j=1, maxneighs*2
                if(update_these(j) == -1) exit
                call update_neighs(m,update_these(j),cutoff)
            enddo
            !call correct_neighs(m,cutoff)
            do jj=1, nneighs(w)
                do j=1, maxneighs*2
                    if(update_these(j) == neighs(w,jj)) exit
                    if(update_these(j) == -1) then
                        update_these(j) = neighs(w,jj)
                        call update_neighs(m,update_these(j),cutoff)
                        exit
                    endif
                enddo
            enddo
            do jj=1, nneighs(w)
                do jjj=1, nneighs(neighs(w,jj))
                    do j=1, maxneighs*2
                        if(update_these(j) == neighs(neighs(w,jj),jjj) ) exit
                        if(update_these(j) == -1) then
                            update_these(j) = neighs(neighs(w,jj),jjj)
                            call update_neighs(m,update_these(j),cutoff)
                            exit
                        endif
                    enddo
                enddo
            enddo
            do jj=1,maxneighs*2
                if(update_these(jj) == -1) exit
                call vp_atom(m,update_these(jj),cutoff)
            enddo
            call group_indexes
            !write(*,*) "New VP index of moved atom:"
            !call print_index(w)
            !write(*,*) neighs(w,:)
            !write(*,*) 'CN:',w,sum(indexes(w,:)),nneighs(w)

            chi2_no_energy = chi_square(alpha,numtypes, numtypes_sim)

            chi2_new = chi2_no_energy + te2
            del_chi = chi2_new - chi2_old
            call mpi_bcast(del_chi, 1, mpi_double, 0, mpi_comm_world, mpierr)

            if(myid .eq. 0) write(*,*) "Energy = ", te2
            if(myid .eq. 0) write(*,*) "LSqF V(k) = ", chi2_no_energy
            if(myid .eq. 0) write(*,*) "cf_old = ", chi2_old
            if(myid .eq. 0) write(*,*) "cf_new = ", chi2_new
            if(myid .eq. 0) write(*,*) "Del-cf = ", del_chi
            if(myid .eq. 0) write(*,*) "Others = ", numtypes_sim(different_types+1)
            if(myid .eq. 0) write(*,*) "Bad = ", numtypes_sim(different_types+2)

            ! Test if the move should be accepted or rejected based on del_chi
            if(del_chi <0.0)then
                ! Accept the move
#ifndef USE_LMP
                e1 = e2 ! eam
#endif
                chi2_old = chi2_new
                accepted = .true.
                if(myid.eq.0) write(*,*) "MC move accepted outright."
            else
                ! Based on the random number above, even if del_chi is negative, decide
                ! whether to move or not (statistically).
                if(log(1.-randnum)<-del_chi*beta)then
                    ! Accept move
#ifndef USE_LMP
                    e1 = e2 ! eam
#endif
                    chi2_old = chi2_new
                    accepted = .true.
                    if(myid.eq.0) write(*,*) "MC move accepted due to probability. del_chi*beta = ", del_chi*beta
                else
                    ! Reject move
                    accepted = .false.
                    if(myid.eq.0) write(*,*) "MC move rejected."
                    call reject_position(m, w, xx_cur, yy_cur, zz_cur)
                    call hutch_move_atom(m,w,xx_cur, yy_cur, zz_cur)  !update hutches.
                    update_these = -1
                    update_these(1) = w
                    do jj=1, nneighs(w)
                        do j=1, maxneighs*2
                            if(update_these(j) == neighs(w,jj)) exit
                            if(update_these(j) == -1) then
                                update_these(j) = neighs(w,jj)
                                exit
                            endif
                        enddo
                    enddo
                    do jj=1, nneighs(w)
                        do jjj=1, nneighs(neighs(w,jj))
                            do j=1, maxneighs*2
                                if(update_these(j) == neighs(neighs(w,jj),jjj) ) exit
                                if(update_these(j) == -1) then
                                    update_these(j) = neighs(neighs(w,jj),jjj)
                                    exit
                                endif
                            enddo
                        enddo
                    enddo
                    do j=1, maxneighs*2
                        if(update_these(j) == -1) exit
                        call update_neighs(m,update_these(j),cutoff)
                    enddo
                    !call correct_neighs(m,cutoff)
                    do jj=1, nneighs(w)
                        do j=1, maxneighs*2
                            if(update_these(j) == neighs(w,jj)) exit
                            if(update_these(j) == -1) then
                                update_these(j) = neighs(w,jj)
                                call update_neighs(m,update_these(j),cutoff)
                                exit
                            endif
                        enddo
                    enddo
                    do jj=1, nneighs(w)
                        do jjj=1, nneighs(neighs(w,jj))
                            do j=1, maxneighs*2
                                if(update_these(j) == neighs(neighs(w,jj),jjj) ) exit
                                if(update_these(j) == -1) then
                                    update_these(j) = neighs(neighs(w,jj),jjj)
                                    call update_neighs(m,update_these(j),cutoff)
                                    exit
                                endif
                            enddo
                        enddo
                    enddo
                    do jj=1,maxneighs*2
                        if(update_these(jj) == -1) exit
                        call vp_atom(m,update_these(jj),cutoff)
                    enddo
                    call group_indexes
                    !write(*,*) "Reverted VP index of moved atom:"
                    !call print_index(w)
                    !write(*,*) 'CN:',w,sum(indexes(w,:)),nneighs(w)
#ifndef USE_LMP
                    e2 = e1 ! eam
#else
                    write(lmp_cmd_str, "(A9, I4, A3, F, A3, F, A3, F)") "set atom ", w, " x ", xx_cur, " y ", yy_cur, " z ", zz_cur
                    call lammps_command(lmp, trim(lmp_cmd_str))
#endif
                endif
            endif
            endif ! if(accepted) from above

            if(myid .eq. 0) then
            if(accepted) then
                acceptance_array(mod(i,100)+1) = 1
            else
                acceptance_array(mod(i,100)+1) = 0
            endif
            if(i .ge. 100) avg_acceptance = sum(acceptance_array)/100.0
            ! Writing to 0 is stderr
            if(i .ge. 100 .and. avg_acceptance .le. 0.05 .and. mod(i,100) .eq. 0) write(0,*) "WARNING!  Acceptance rate is low:", avg_acceptance
            endif

            ! Periodically save data.
            if(myid .eq. 0) then
            if(mod(i,2000)==0)then
                ! Write to vk_update ERROR HERE - if the most recent move was
                ! rejected then this will print the incorrect vk. TODO
                !write(vku_fn, "(A9)") "vk_update"
                !write(step_str,*) i
                !vku_fn = trim(trim(trim(trim(vku_fn)//jobID)//"_")//step_str)//".txt"
                !open(32,file=trim(vku_fn),form='formatted',status='unknown')
                !    do j=1, nk
                !        write(32,*)k(j),vk(j)
                !    enddo
                !close(32)
                ! Write to model_update
                ! This takes a bit of time.
                write(output_model_fn, "(A12)") "model_update"
                write(step_str,*) i
                output_model_fn = trim(trim(trim(trim(output_model_fn)//jobID)//"_")//step_str)//".xyz"
                open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                    write(33,*)"updated model"
                    write(33,*)m%lx,m%ly,m%lz
                    do j=1,m%natoms
                        write(33,*)m%znum%ind(j), m%xx%ind(j), m%yy%ind(j), m%zz%ind(j)
                    enddo
                    write(33,*)"-1"
                    !do j=1,m%natoms
                    !    write(33,'(9I3)') j, indexes(j,:)
                    !enddo
                close(33)
                !do j=myid+1,nrot,numprocs
                !    write(myid_str,*) j
                !    !output_model_fn = trim(trim(trim(trim(trim(trim(output_model_fn)//jobID)//"_")//step_str)//"_")//myid_str)//".xyz"
                !    output_model_fn = trim(trim("mrot_model")//trim(myid_str))
                !    !write(*,*) trim(output_model_fn)
                !    open(33,file=trim(output_model_fn),form='formatted',status='unknown')
                !        write(33,*)"mrot model. rot=", j, "step=",i
                !        write(33,*)m%lx,m%ly,m%lz
                !        do ii=1,mrot(j)%natoms
                !            write(33,*)mrot(j)%znum%ind(ii), mrot(j)%xx%ind(ii), mrot(j)%yy%ind(ii), mrot(j)%zz%ind(ii)
                !        enddo
                !        write(33,*)"-1"
                !    close(33)
                !enddo
            endif
            if(mod(i,1)==0)then
                !if(accepted) then
                    ! Write chi2 info
                    open(36,file=trim(chi_squared_file),form='formatted',status='unknown',access='append')
                        write(36,*) i, chi2_no_energy, te2
                    close(36)
                !endif
            endif
#ifdef TIMING
            if(mod(i,1)==0)then
                ! Write to time_elapsed
                open(35,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
                    t1 = omp_get_wtime()
                    write (35,*) i, t1-t0, (t1-t0)/i, t1-t2
                close(35)
            endif
#endif
            if(mod(i,100)==0 .and. i .ge. 100)then
                ! Write to acceptance rate
                open(40,file=trim(acceptance_rate_fn),form='formatted',status='unknown',access='append')
                    write(40,*) i, avg_acceptance
                close(40)
            endif
            endif ! myid == 0

            ! Every 400,000 steps lower the temp, max_move, and reset beta.
            if(mod(i,400000)==0)then
                temperature = temperature * sqrt(0.7)
                if(myid.eq.0) write(*,*) "Lowering temp to", temperature, "at step", i
                max_move = max_move * sqrt(0.94)
                beta=1./((8.6171e-05)*temperature)
            endif
            
            i=i+1

        enddo !RMC do loop
        write(*,*) "Monte Carlo Finished!"

        ! The rmc loop finished. Write final data.
        if(myid.eq.0)then
            ! Write final model
            open(unit=55,file=trim(final_model_fn),form='formatted',status='unknown')
            write(55,*)"final model"
            write(55,*)m%lx,m%ly,m%lz
            do i=1,m%natoms
                write(55,*)m%znum%ind(i), m%xx%ind(i), m%yy%ind(i), m%zz%ind(i)
            enddo
            write(55,*)"-1"; close(55)
            ! Write final energy.
            !open(56,file=trim(energy_fn),form='formatted', status='unknown',access='append')
            !write(56,*) i, te2
            !close(56)
#ifdef TIMING
            ! Write final time spent.
            open(57,file=trim(time_elapsed),form='formatted',status='unknown',access='append')
            t1 = omp_get_wtime()
            write (57,*) i, t1-t0
            write(57,*) "Finshed.", numprocs, "processors."
            close(57)
#endif
        endif
    endif ! Use RMC

#ifdef USE_LMP
    if(myid.eq.0) then
    call lammps_close (lmp)
    endif
#endif
    call mpi_finalize(mpierr)

end program rmc
