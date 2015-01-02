

module vor_mod
    use model_mod
    ! Move stuff up above into here
    ! Include bad, indexes, neighs
    implicit none
    integer, parameter :: n0 = 50000 ! Max number of atoms?
    integer, parameter :: maxcan = 75 !75
    integer, parameter :: maxver = 100 !100
    integer, parameter :: maxepf = 20  !20
    integer, parameter :: maxneighs = 50
    real, parameter :: tol = 0.03
    real, save :: p(maxcan,4),v(maxver,3) ! vertices
    integer, allocatable, dimension(:), save :: bad, nneighs
    integer, allocatable, dimension(:,:), save :: indexes
    integer, allocatable, dimension(:,:), save :: neighs
    integer, allocatable, dimension(:,:), save :: types
    integer, allocatable, dimension(:), save :: numtypes, numtypes_sim
    integer, save :: mtag(maxcan) ! tag of atom after sorting
    integer, save :: nloop(maxepf,maxcan),nepf(maxcan) ! check connection
    integer, save :: mvijk(maxver,3)
    integer, save :: nc, nv, ne, nf ! # of candidates,vertices,edges,faces
    integer, save :: natoms, nlines
    integer, save :: different_types

    contains
    subroutine vor_init(m, rcut, paramfile)
        ! Intializes data that will be saves for each atom
        ! Intializes neighbors for each atom, max neighbors is set here
        implicit none
        type(model), intent(in) :: m ! Model
        real, intent(in) :: rcut
        integer :: i, j
        integer :: nlist, istat
        integer, pointer, dimension(:) :: atoms
        real :: x2, y2, z2, rr2, rcut2
        integer :: numneighs
        character (len=256), intent(in) :: paramfile

        natoms = m%natoms
        rcut2 = rcut**2
        allocate(nneighs(m%natoms))
        nneighs = 0
        allocate(bad(m%natoms))
        bad = 0
        allocate(indexes(m%natoms,8))
        indexes = 0

        allocate(neighs(m%natoms,maxneighs))
        ! Find neighbors
        do i=1, m%natoms
            call hutch_list_3D(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), rcut, atoms, istat, nlist)
            numneighs = 0
            do j=1, nlist-1
                if(i .ne. atoms(j)) then
                    x2 = m%xx%ind(atoms(j)) - m%xx%ind(i)
                    y2 = m%yy%ind(atoms(j)) - m%yy%ind(i)
                    z2 = m%zz%ind(atoms(j)) - m%zz%ind(i)
                    x2 = x2 - m%lx*anint(x2/m%lx)
                    y2 = y2 - m%ly*anint(y2/m%ly)
                    z2 = z2 - m%lz*anint(z2/m%lz)
                    rr2 = x2**2 + y2**2 + z2**2
                    if( rr2 < rcut2 ) then
                        numneighs = numneighs + 1
                        neighs(i,numneighs) = atoms(j)
                    endif
                endif
            enddo
            if(associated(atoms)) deallocate(atoms)
            do j=numneighs+1, maxneighs
                neighs(i,j) = 0
            enddo
            nneighs(i) = numneighs
        enddo

        !call group_indexes_from_paramfile('/home/jjmaldonis/development/descriptors/vprmc/src/temp.txt')
        call group_indexes_from_paramfile(paramfile)
    end subroutine vor_init

    subroutine update_neighs(m, atom, rcut)
        ! Updates the neighbors of one atom
        implicit none
        type(model), intent(in) :: m ! Model
        integer, intent(in) :: atom
        real, intent(in) :: rcut
        integer :: nlist, istat
        integer :: j
        integer, pointer, dimension(:) :: atoms
        real :: x2, y2, z2, rr2, rcut2
        integer :: numneighs
        rcut2 = rcut**2
        call hutch_list_3D(m, m%xx%ind(atom), m%yy%ind(atom), m%zz%ind(atom), rcut, atoms, istat, nlist)
        numneighs = 0
        do j=1, nlist-1
            if(atom .ne. atoms(j)) then
                x2 = m%xx%ind(atoms(j)) - m%xx%ind(atom)
                y2 = m%yy%ind(atoms(j)) - m%yy%ind(atom)
                z2 = m%zz%ind(atoms(j)) - m%zz%ind(atom)
                x2 = x2 - m%lx*anint(x2/m%lx)
                y2 = y2 - m%ly*anint(y2/m%ly)
                z2 = z2 - m%lz*anint(z2/m%lz)
                rr2 = x2**2 + y2**2 + z2**2
                !write(*,*) atoms(j), x2,y2,z2, sqrt(rr2), rcut
                if( rr2 < rcut2 ) then
                    numneighs = numneighs + 1
                    neighs(atom,numneighs) = atoms(j)
                endif
            endif
        enddo
        if(associated(atoms)) deallocate(atoms)
        do j=numneighs+1, maxneighs
            neighs(atom,j) = 0
        enddo
        nneighs(atom) = numneighs
        !write(*,*) "New neighs for atom",atom,nneighs(atom)
        !write(*,*) neighs(atom,1:numneighs)
        !write(*,*) numneighs
        !write(*,*) neighs(atom,:)
    end subroutine update_neighs

    subroutine vp_atom(m,atom,rcut)
        ! Calculates the VP for one atom
        implicit none
        type(model), intent(in) :: m ! Model
        integer, intent(in) :: atom
        real, intent(in) :: rcut ! Cutoff
        integer :: j
        integer :: ic
        real :: rcutsq
        real :: rxij, ryij, rzij, rijsq

        bad(atom) = 0
        ic = 0
        rcutsq=rcut**2
        do j=1, m%natoms
            if(j .ne. atom) then
                rxij = (m%xx%ind(j)+m%lx/2.0)/m%lx - (m%xx%ind(atom)+m%lx/2.0)/m%lx
                ryij = (m%yy%ind(j)+m%ly/2.0)/m%ly - (m%yy%ind(atom)+m%ly/2.0)/m%ly
                rzij = (m%zz%ind(j)+m%lz/2.0)/m%lz - (m%zz%ind(atom)+m%lz/2.0)/m%lz
                rxij = rxij - anint ( rxij )
                ryij = ryij - anint ( ryij )
                rzij = rzij - anint ( rzij )
                rxij=rxij*m%lx
                ryij=ryij*m%ly
                rzij=rzij*m%lz
                rijsq = rxij**2 + ryij**2 + rzij**2

                if ( rijsq .lt. rcutsq ) then
                    ic = ic + 1
                    if(ic.gt.maxcan)then
                        write(*,*)ic, maxcan
                        stop 'too many candidates'
                    endif
                    p(ic,1) = rxij
                    p(ic,2) = ryij
                    p(ic,3) = rzij
                    p(ic,4) = rijsq
                    mtag(ic) = j
                endif
            endif
        enddo

        nc = ic
        call sort_pmtag
        call work(atom)
        if(bad(atom)) then
            do j=1,8
                indexes(atom,j) = 0
            enddo
            write(*,*) "Atom", atom, "is bad!"
        else
            call vpindex(atom)
        endif
    end subroutine vp_atom

    subroutine work(atom)
        implicit none
        integer, intent(in) :: atom
        integer :: i, j, k, iv, l
        real :: ai, bi, ci, di, aj, bj, cj, dj, ab, bc, ca, da, db, dc, ak, bk, ck, dk
        real :: det, detinv, vxijk, vyijk, vzijk
        logical :: ok

        if(nc.lt.4) then
            write(*,*) 'less than 4 points given to work', nc
            bad(atom) = 1
            return
            !stop
        endif

        iv = 0
        do i=1, nc-2
            ai = p(i,1)
            bi = p(i,2)
            ci = p(i,3)
            di = -p(i,4)
            do j=i+1, nc-1
                aj =  p(j,1)
                bj =  p(j,2)
                cj =  p(j,3)
                dj = -p(j,4)
                ab = ai * bj - aj * bi
                bc = bi * cj - bj * ci
                ca = ci * aj - cj * ai
                da = di * aj - dj * ai
                db = di * bj - dj * bi
                dc = di * cj - dj * ci
                do k=j+1, nc
                    ak =  p(k,1)
                    bk =  p(k,2)
                    ck =  p(k,3)
                    dk = -p(k,4)
                    det = ak * bc + bk * ca + ck * ab
                    if ( abs ( det ) .gt. tol ) then
                        detinv = 1.0 / det
                        vxijk = ( - dk * bc + bk * dc - ck * db ) * detinv
                        vyijk = ( - ak * dc - dk * ca + ck * da ) * detinv
                        vzijk = (   ak * db - bk * da - dk * ab ) * detinv
                        ok = .true.
                        l=1
                        do while( ok .and. (l .le. nc) )
                            if((l.ne.i).and.(l.ne.j).and.(l.ne.k)) then
                                ok=((p(l,1)*vxijk+p(l,2)*vyijk+p(l,3)*vzijk).le.p(l,4))
                            endif
                            l=l+1
                        enddo
                        if(ok) then
                            iv = iv + 1
                            if(iv .gt. maxver) stop 'too many vertices'
                            mvijk(iv,1)  = i
                            mvijk(iv,2)  = j
                            mvijk(iv,3)  = k
                            v(iv,1) = 0.5 * vxijk
                            v(iv,2) = 0.5 * vyijk
                            v(iv,3) = 0.5 * vzijk
                            !write(*,*) "DEBUG2",iv,vxijk,vyijk,vzijk
                        endif
                    endif
                enddo
            enddo
        enddo

        nv = iv
        if(nv.lt.4) then
            bad(atom) = 1
            return
            !stop 'less than 4 vertices found in work'
        endif

        ! These are arrays
        nepf = 0
        nloop = 0
        !write(*,*) "Number of vertices = ", nv
        !call sleep(1)
        do iv=1, nv ! mvijk is set in the loop above.
            nepf(mvijk(iv,1)) = nepf(mvijk(iv,1)) + 1
            nepf(mvijk(iv,2)) = nepf(mvijk(iv,2)) + 1
            nepf(mvijk(iv,3)) = nepf(mvijk(iv,3)) + 1
            !write(*,*) nepf(mvijk(iv,1)),nepf(mvijk(iv,2)),nepf(mvijk(iv,3)) !TODO
            if(nepf(mvijk(iv,1)).gt.maxepf) stop 'epf>maxepf'
            if(nepf(mvijk(iv,2)).gt.maxepf) stop 'epf>maxepf'
            if(nepf(mvijk(iv,3)).gt.maxepf) stop 'epf>maxepf'
            nloop(nepf(mvijk(iv,1)),mvijk(iv,1)) = iv
            nloop(nepf(mvijk(iv,2)),mvijk(iv,2)) = iv
            nloop(nepf(mvijk(iv,3)),mvijk(iv,3)) = iv
            !write(*,*)
            !nepf(mvijk(iv,1)),mvijk(iv,1),nloop(nepf(mvijk(iv,1)),mvijk(iv,1))
        enddo
        !write(*,*) nloop

        nf = 0
        ne = 0
        do i=1, nc
            if(nepf(i) .gt. 0) nf = nf + 1
            ne = ne + nepf(i)
        enddo
        !if(mod(ne,2).ne.0)then
        !    write(*,*)"ne=", ne
        !    do iv=1,nv
        !        write(*,*)"mvijk(iv,:)=", mvijk(iv,:)
        !    enddo
        !endif
        ne = ne/2
        if((nv-ne+nf).ne.2)then
            bad(atom) = 1
        endif
    end subroutine work


    subroutine sort_pmtag!(nc, p, mtag)
        implicit none
        integer :: itag, i, j
        real :: pi1, pi2, pi3, pi4 ! Temp vars for sorting purposes.
        do i=1,nc-1
            do j=i+1,nc
                if(p(i,4).gt.p(j,4))then
                    pi1 = p(i,1)
                    pi2 = p(i,2)
                    pi3 = p(i,3)
                    pi4 = p(i,4)
                    itag = mtag(i)
                    p(i,1) = p(j,1)
                    p(i,2) = p(j,2)
                    p(i,3) = p(j,3)
                    p(i,4) = p(j,4)
                    mtag(i) = mtag(j)
                    p(j,1) = pi1
                    p(j,2) = pi2
                    p(j,3) = pi3
                    p(j,4) = pi4
                    mtag(j) = itag
                end if
            enddo
        enddo
    end subroutine sort_pmtag

    subroutine vpindex(atom)
        ! Converts nepf to an index and saves it into indexes
        implicit none
        integer, intent(in) :: atom
        integer :: i
        !write(*,*) atom
        !write(*,*) nepf(1:nc)
        do i=1,8
            indexes(atom,i) = 0
        enddo
        do i=1, nneighs(atom)
            if(nepf(i) > 0) then
                indexes(atom,nepf(i)-2) = indexes(atom,nepf(i)-2) + 1
            endif
        enddo
    end subroutine vpindex

  SUBROUTINE GetLine(unit, line, stat, iomsg)   
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_EOR   
    !---------------------------------------------------------------------------
    ! Arguments   
    INTEGER, INTENT(IN) :: unit
    CHARACTER(len=256), INTENT(OUT) :: line
    INTEGER, INTENT(OUT) :: stat
    CHARACTER(*), INTENT(OUT) :: iomsg   
    !---------------------------------------------------------------------------
    ! Local variables   
    ! Buffer to read the line (or partial line).
    CHARACTER(256) :: buffer   
    INTEGER :: size           ! Number of characters read from the file.   
    !***************************************************************************   
    line = ''
    DO
      READ (unit, "(A)", ADVANCE='NO', IOSTAT=stat, IOMSG=iomsg, SIZE=size)  &
          buffer
      IF (stat > 0) RETURN      ! Some sort of error.
      line = line // buffer(:size)
      IF (stat < 0) THEN
        IF (stat == IOSTAT_EOR) stat = 0
        RETURN
      END IF
    END DO   
  END SUBROUTINE GetLine

    subroutine group_indexes_from_paramfile(paramfile)
        implicit none
        !character (len=*), intent(in) :: paramfile
        character (len=*) :: paramfile
        integer :: i, j, k, ii
        integer, allocatable, dimension(:,:) :: temp_types
        logical :: found
        integer :: stat
        CHARACTER (len=256) :: line, iomsg

        allocate(temp_types(natoms,8))
        open(12,file=trim(paramfile),form='formatted',status='unknown')

        ii = 0
        do while(.true.)
            call GetLine(12,line,stat,iomsg)
            if(stat < 0) exit
            ii = ii + 1
        enddo
        rewind(12)

        nlines = ii
        do i=1, nlines
            read(12, *) temp_types(i,:)
        enddo
        close(12)

        different_types = 1
        do i=2, nlines! skip the first index, we know its unique
            found = .false.
            do k=1, i-1
                found = .true.
                do j=1, 8
                    if( temp_types(i,j) .ne. temp_types(k,j) ) found = .false.
                enddo
                if(found) exit
            enddo
            if(.not. found) different_types = different_types + 1
        enddo
        write(*,*) "Diff types:", different_types
        allocate(types(different_types,8))
        types = 0
        allocate(numtypes(different_types+2)) ! The first extra spot is for other indexes not in types, the second is for bad atoms
        numtypes = 0
        allocate(numtypes_sim(different_types+2))
        numtypes_sim = 0
        ii = 1
        do j=1,8
            types(ii,j) = temp_types(1,j)
        enddo
        do i=2, nlines! skip the first index, we know its unique
            found = .false.
            do k=1, i-1
                found = .true.
                do j=1, 8
                    if( temp_types(i,j) .ne. temp_types(k,j) ) found = .false.
                enddo
                if(found) exit
            enddo
            if(.not. found) then
                ii = ii + 1
                do j=1, 8
                    types(ii,j) = temp_types(i,j)
                enddo
            endif
        enddo

        do i=1, natoms
            if(bad(i)) then
                numtypes(different_types+2) = numtypes(different_types+2) + 1
            else
                do k=1, different_types
                    found = .true.
                    do j=1, 8
                        if( temp_types(i,j) .ne. types(k,j) ) found = .false.
                    enddo
                    if(found) exit
                enddo
                if(found) then
                    numtypes(k) = numtypes(k) + 1
                else
                    numtypes(different_types+1) = numtypes(different_types+1) + 1
                endif
            endif
        enddo

        !do k=1,different_types
        !    write(*,'(I4,A3,8I3)') numtypes(k), '   ', types(k,:)
        !enddo
        !write(*,*) sum(numtypes)
        !write(*,*) numtypes

    end subroutine group_indexes_from_paramfile

    subroutine group_indexes
        integer :: i, j, k
        logical :: found
        do k=1, different_types+1
            numtypes_sim(k) = 0
        enddo
        do i=1, natoms
            if(bad(i)) then
                numtypes(different_types+2) = numtypes(different_types+2) + 1
            else
                do k=1, different_types
                    found = .true.
                    do j=1, 8
                        if( indexes(i,j) .ne. types(k,j) ) found = .false.
                    enddo
                    if(found) exit
                enddo
                if(found) then
                    numtypes_sim(k) = numtypes_sim(k) + 1
                else
                    numtypes_sim(different_types+1) = numtypes_sim(different_types+1) + 1
                endif
            endif
        enddo
        !do k=1,different_types
        !    write(*,'(I4,A3,8I3)') numtypes_sim(k), '   ', types(k,:)
        !enddo
        !write(*,*) sum(numtypes_sim)
        !write(*,*) numtypes_sim(different_types+1)
        !write(*,*) numtypes_sim
    end subroutine group_indexes
    
    subroutine print_index(atom)
        integer, intent(in) :: atom
        integer :: j
        character (len=32) :: cc, s

        s = ''
        do j=1,8
            write(cc,'(I2)') indexes(atom,j)
            s = trim(s)//' '//trim(cc)
        enddo
        write(*,*) trim(s)
        !write(*,*) neighs(atom,:)
    end subroutine print_index

    subroutine check_neighs(m, rcut)
        ! Goes through neighs and nneighs and makes sure the values are correct
        implicit none
        type(model), intent(in) :: m ! Model
        real, intent(in) :: rcut
        integer :: i, j
        integer :: nlist, istat
        integer, pointer, dimension(:) :: atoms
        real :: x2, y2, z2, rr2, rcut2
        integer :: numneighs

        natoms = m%natoms
        rcut2 = rcut**2

        ! Find neighbors
        do i=1, m%natoms
            call hutch_list_3D(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), rcut, atoms, istat, nlist)
            numneighs = 0
            do j=1, nlist-1
                if(i .ne. atoms(j)) then
                    x2 = m%xx%ind(atoms(j)) - m%xx%ind(i)
                    y2 = m%yy%ind(atoms(j)) - m%yy%ind(i)
                    z2 = m%zz%ind(atoms(j)) - m%zz%ind(i)
                    x2 = x2 - m%lx*anint(x2/m%lx)
                    y2 = y2 - m%ly*anint(y2/m%ly)
                    z2 = z2 - m%lz*anint(z2/m%lz)
                    rr2 = x2**2 + y2**2 + z2**2
                    if( rr2 < rcut2 ) then
                        numneighs = numneighs + 1
                        if(neighs(i,numneighs) .ne. atoms(j)) write(*,*) "HOUSTON! neighs is wrong!",i, neighs(i,numneighs), atoms(j)
                    endif
                endif
            enddo
            if(associated(atoms)) deallocate(atoms)
            do j=numneighs+1, maxneighs
                if(neighs(i,j) .ne. 0) write(*,*) "HOUSTON! No neigh should be here!",i, neighs(i,j)
            enddo
            if(nneighs(i) .ne. numneighs) write(*,*) "HOUSTON! nneighs not equal!",i, nneighs(i), numneighs
        enddo
    end subroutine check_neighs


    subroutine correct_neighs(m, rcut)
        ! Goes through neighs and nneighs and makes sure the values are correct
        implicit none
        type(model), intent(in) :: m ! Model
        real, intent(in) :: rcut
        integer :: i, j
        integer :: nlist, istat
        integer, pointer, dimension(:) :: atoms
        real :: x2, y2, z2, rr2, rcut2
        integer :: numneighs

        natoms = m%natoms
        rcut2 = rcut**2

        ! Find neighbors
        do i=1, m%natoms
            call hutch_list_3D(m, m%xx%ind(i), m%yy%ind(i), m%zz%ind(i), rcut, atoms, istat, nlist)
            numneighs = 0
            do j=1, nlist-1
                if(i .ne. atoms(j)) then
                    x2 = m%xx%ind(atoms(j)) - m%xx%ind(i)
                    y2 = m%yy%ind(atoms(j)) - m%yy%ind(i)
                    z2 = m%zz%ind(atoms(j)) - m%zz%ind(i)
                    x2 = x2 - m%lx*anint(x2/m%lx)
                    y2 = y2 - m%ly*anint(y2/m%ly)
                    z2 = z2 - m%lz*anint(z2/m%lz)
                    rr2 = x2**2 + y2**2 + z2**2
                    if( rr2 < rcut2 ) then
                        numneighs = numneighs + 1
                        neighs(i,numneighs) = atoms(j)
                    endif
                endif
            enddo
            if(associated(atoms)) deallocate(atoms)
            do j=numneighs+1, maxneighs
                neighs(i,j) = 0
            enddo
            nneighs(i) = numneighs
        enddo
    end subroutine correct_neighs


end module vor_mod








!program main
!    use vor_mod
!    use model_mod
!    implicit none
!    type(model) :: m
!    character (len=256) :: model_filename, c
!    integer :: istat, length, i, j
!    real :: cutoff
!    character (len=32) :: cc, s
!
!    call get_command_argument(1, c, length, istat)
!    if (istat == 0) then
!        model_filename = trim(c)
!    else
!        error stop "istat for model_filename get_command_arg was nonzero"
!    endif
!    model_filename = trim(model_filename)
!    call read_model(model_filename, m, istat)
!    call check_model(m, istat)
!    call recenter_model(0.0, 0.0, 0.0, m)
!
!    cutoff = 3.6
!    call vor_init(m,cutoff)
!    do i=1, m%natoms
!        call vp_atom(m, i, cutoff)
!        s = ''
!        do j=1,8
!            write(cc,'(I2)') indexes(i,j)
!            s = trim(s)//' '//trim(cc)
!        enddo
!        !write(*,*) i, trim(s)
!    enddo
!end program  main




