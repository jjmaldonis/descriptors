MODULE rmc_functions

  use model_mod
  use vor_mod
  implicit none

  interface
     function ran2(idum)
       real :: ran2
       integer :: idum
     end function ran2
  end interface

CONTAINS

    function chi_square(alpha,vpexp, vpsim)
        integer, allocatable, dimension(:) :: vpexp, vpsim
        real, intent(in) :: alpha
        integer :: i, valexp, valsim, others
        double precision :: chi_square
        chi_square=0.0
        !vpexp = vpexp*8

        others = vpsim(different_types+1)
        do i=1,different_types
            valsim = vpsim(i)
            valexp = vpexp(i)
            if(valsim > valexp) then
                others = others + (valsim-valexp)
                valsim = valexp
            endif
            chi_square = chi_square + (valsim - valexp)**2/sqrt(float(vpexp(i)))
            !write(*,*) i, chi_square
        enddo
        chi_square = chi_square + (others-vpexp(different_types+1))**2/sqrt(float(vpexp(different_types+1))) ! This is for "others", ignore the stuff below

        !chi_square = chi_square/(different_types+1) ! Commenting this because I
        ! don't want to penalize for not having the exact amount of "others". It
        ! is okay for the "others" be in a different category.

        chi_square = chi_square + vpsim(different_types+2)/(different_types+1)*3.0 ! Weight the "bad" atoms as 3 times as bad to have
        chi_square = chi_square * alpha
    end function chi_square

    function check_cutoffs(m,cutoff_r,moved_atom)
        logical check_cutoffs
        real,dimension(:,:) ::cutoff_r
        integer moved_atom
        type(model) m
        integer, dimension(:), pointer :: atoms
        integer  nlist
        integer istat
        real radius, temp_x, temp_y, temp_z
        integer i,j
        integer num1, num2
        real dist_pair

        !find the maximum cut-off distance
        radius=maxval(maxval(cutoff_r, 1),1)

        ! Find all atoms within radius radius of moved_atom and put them 
        ! in the list atoms. Also sets nlist == size(atoms)+1.
        call hutch_list_3D(m, m%xx%ind(moved_atom),m%yy%ind(moved_atom),m%zz%ind(moved_atom), radius, atoms, istat, nlist)
        if (istat .eq. 1) then
            print *, 'memory allocation fails!'
            return
        else if(istat .eq. -1) then
            print *, 'no atom is found', radius
            check_cutoffs = .true. 
            return
        endif

        !begin to calculate pair distance
        !note, there are (nlist-1) atoms in total in 'atoms'

        !first, determine the type of moved_atom 
        do i=1, m%nelements
            if(m%znum%ind(moved_atom) .eq. m%atom_type(i)) then
                num1 = i
                exit
            endif
        enddo
        
        do i=1, (nlist-1)
            !do not count the pair distance to itself
            ! The below if statement should be irrelevant. moved_atom isn't put
            ! in atoms by hutch_list_3d as far as i know.
            check_cutoffs = .true.  !jwh 032409
            if (atoms(i) .ne. moved_atom) then
                do j=1, m%nelements
                    if(m%atom_type(j) .eq. m%znum%ind(atoms(i))) then
                        num2 = j
                        exit
                    endif
                enddo !j=1, m%nelements
              
                !endif  !032409 - jwh
               
                !calculate the atomic distance
                !compare with cutoff_r
                temp_x = m%xx%ind(moved_atom) - m%xx%ind(atoms(i))  !pbc added - jwh 04/14/2009
                temp_y = m%yy%ind(moved_atom) - m%yy%ind(atoms(i))
                temp_z = m%zz%ind(moved_atom) - m%zz%ind(atoms(i))
                temp_x = temp_x - m%lx*anint(temp_x/m%lx)
                temp_y = temp_y - m%ly*anint(temp_y/m%ly)
                temp_z = temp_z - m%lz*anint(temp_z/m%lz)
                  
                dist_pair = temp_x**2 + temp_y**2 + temp_z**2
                dist_pair = sqrt(dist_pair)
                   
                if (dist_pair  .lt. cutoff_r(num1, num2)) then
                    !write(*,*)"DEBUG", dist_pair  , cutoff_r(num1, num2) !debug - jwh 032409
                    check_cutoffs=.false.
                    exit
                endif
            endif !032409 - jwh
        enddo !i=1, (nlist-1)

        !if (associated(atoms)) deallocate(atoms) ! pmv 4/17/09
        if (nlist .gt. 1) deallocate(atoms) ! fy 4/17/09
    end function check_cutoffs

    subroutine random_move(m,w,xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new, alpha)
        use mpi  
        type(model), intent(inout) :: m
        integer iseed, w
        real alpha, aa, bb, cc
        real, intent(out) :: xx_cur,yy_cur,zz_cur,xx_new,yy_new,zz_new    !Cur and new positions of atoms
        real :: rand1, rand2, rand3, rand4


        iseed = 791315

        rand1 = ran2(iseed)
        rand2 = ran2(iseed)
        rand3 = ran2(iseed)
        rand4 = ran2(iseed)

        w = int(m%natoms*rand1)+1

        !write(*,*)myid, rand1, rand2, rand3, rand4

        xx_cur = m%xx%ind(w)            !Cur positions of the atom before random move
        yy_cur = m%yy%ind(w)
        zz_cur = m%zz%ind(w)
        
        aa = alpha*(rand2 - 0.5)
        bb = alpha*(rand3 - 0.5)
        cc = alpha*(rand4 - 0.5)
        
        m%xx%ind(w) = m%xx%ind(w) + aa 
        m%yy%ind(w) = m%yy%ind(w) + bb 
        m%zz%ind(w) = m%zz%ind(w) + cc
        
        if(m%xx%ind(w)>m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)-m%lx       !pbc 
        if(m%yy%ind(w)>m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)-m%ly
        if(m%zz%ind(w)>m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)-m%lz
        if(m%xx%ind(w)<-m%lx*0.5) m%xx%ind(w)=m%xx%ind(w)+m%lx
        if(m%yy%ind(w)<-m%ly*0.5) m%yy%ind(w)=m%yy%ind(w)+m%ly
        if(m%zz%ind(w)<-m%lz*0.5) m%zz%ind(w)=m%zz%ind(w)+m%lz
        
        xx_new=m%xx%ind(w)              !new positions of the atom after random move
        yy_new=m%yy%ind(w)
        zz_new=m%zz%ind(w)
    end subroutine random_move


END MODULE rmc_functions
