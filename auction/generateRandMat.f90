PROGRAM main

	implicit none

	include 'mpif.h'
	include 'omp_lib.h'

	integer :: mat_size
	integer :: rank, nProcs, error
	integer, allocatable :: ia(:), ja(:)
	double precision, allocatable :: a(:), tmp_a(:)
	integer, allocatable :: rBlockBndry(:)
	double precision :: option(13)
	integer, allocatable :: match_perm(:)
	double precision, allocatable :: match_scaling(:)
	integer :: i, j
	character*20 :: ineqns
	integer :: seed(8)
	integer :: nCols, tot_neqns, neqns, nnz
	integer*8 :: tot_nnz
	integer :: parS 
	double precision :: sumVal

	integer :: mc64_, auction_, greedy_
	character*1 :: mmc64, mauction, mgreedy
	integer(kind=1) :: iargc

	neqns = 0; nnz = 0; tot_neqns = 0; tot_nnz = 0
	i = 0; j = 0; error = 0
	mc64_ = 0; auction_ = 1; greedy_ = 0

	parS = 1

	! MPI Settings
	call MPI_INIT (error)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, nProcs, error)
	call MPI_COMM_RANK (MPI_COMM_WORLD, rank, error)
	call MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN, error)

	if(rank .eq. 0) then
		call getarg(1, ineqns)
		ineqns = trim(ineqns)
 		read(ineqns, '(i10)') neqns
		write(*,*) 'matrix local neqns', neqns
		if(parS .eq. 0) then
			nCols = neqns * nProcs
			neqns = neqns / nProcs
		else
			nCols = neqns
			neqns = neqns / nProcs
		end if
		
		if(iargc() > 1) then
			call getarg (2, mmc64)
			read(mmc64, '(i10)') mc64_
		end if
		if(iargc() > 2) then
			call getarg (3, mauction)
			read(mauction, '(i10)') auction_
		end if
		if(iargc() > 3) then
			call getarg (4, mgreedy)
			read(mgreedy, '(i10)') greedy_
		end if

	end if
	
	call MPI_BCAST(nCols, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
	call MPI_BCAST(neqns, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
	call MPI_BCAST(auction_, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)

	nnz = neqns * nCols
	write(*,*) 'nrows, ncols, nnz, rank', neqns, nCols, nnz, rank, neqns*nCols

	allocate(ia(neqns+1))
	allocate(ja(nnz))
	allocate(a(nnz))
	!!allocate(tmp_a(nnz*nProcs))

	ia = 0; ja = 0; a = 0.0d0

	seed(1) = 200 / (rank+1)
	seed(2) = 2000 / (rank+1)

	call RANDOM_SEED(put=seed)
	!!call RANDOM_SEED

	!!if(rank .eq. 0) then
	!!	do i = 1, nnz*nProcs
	!!		call random_number(tmp_a(i))
	!!	end do
	!!end if

	do i = 1, nnz
		call random_number(a(i))
	end do

!!	call MPI_BCAST(tmp_a, nnz*nProcs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)

!!	do i = 1, nnz
!!		a(i) = tmp_a(nnz*rank+i)
!!	end do

	!!deallocate(tmp_a)

	do i = 0, neqns
		ia(i+1) = i*nCols
	end do

	do i = 1, neqns
		do j = 0, nCols-1
			ja((i-1)*nCols+j+1) = j
		end do
	end do
#if 0
    write(10+rank,*) '%%MatrixMarket matrix coordinate real general'
	write(10+rank,*) neqns, neqns, nnz
	do i = 1, neqns
       do j = ia(i)+1, ia(i+1)
          write(10+rank,*) i, ja(j)+1, a(j)
       end do    
    end do
   
	call MPI_Barrier(MPI_COMM_WORLD)
	stop
#endif
#if 0
	write(10+rank,*) neqns
	write(10+rank,*) nnz
	do i = 1, neqns + 1
		write(10+rank,*) ia(i)
	end do
	do i = 1, nnz
		write(10+rank,*) ja(i)
	end do
	do i = 1, nnz
		write(10+rank,*) a(i)
	end do

	call MPI_Barrier(MPI_COMM_WORLD)
	stop
#endif
	tot_neqns = nCols
	tot_nnz = tot_neqns * tot_neqns

	write(*,*) neqns, nnz, tot_neqns, tot_nnz
	!!call MPI_Barrier(MPI_COMM_WORLD)
	!!stop

	if(rank .eq. 0) then

		allocate(match_perm(tot_neqns))
		!! not required
		allocate(match_scaling(1))
		allocate(rBlockBndry(nProcs))

		match_perm = 0; match_scaling = 0
		
		rBlockBndry(1) = 0

		do i = 1, nProcs-1
			rBlockBndry(i+1) = rBlockBndry(i) + neqns
		end do
		
        	if(ia(1) .eq. 1) then
        		!$OMP PARALLEL DO PRIVATE(i)
        		do i = 1, neqns+1
        			ia(i) = ia(i)-1
        		end do
        		!$OMP END PARALLEL DO
        		
        		!$OMP PARALLEL DO PRIVATE(i)
        		do i = 1, nnz
        			ja(i) = ja(i)-1
        		end do
        		!$OMP END PARALLEL DO
        	end if
		
		
	end if

    option(1) = 1 !! info
   	option(2) = 0 !! factor
  	option(3) = 0 !! theta
 	option(4) = 0 !! omega
  	option(5) = 0 !! iter_local
	option(6) = 0 !! iter_all
    option(7) = 1 !! approx_factor
   	option(8) = 1 !! matrix distributed
  	option(9) = 1 !! boundary info
 	option(10) = 1 !! dense
  	option(11) = 0 !! return dual var
	option(12) = 0 !! eps scaling
	option(13) = 0 !! values scaled

	if(auction_ .eq. 1) then
	
	    allocate(tmp_a(nnz))
	    
    	!$OMP PARALLEL DO PRIVATE(i)
	    do i = 1, nnz
	    	tmp_a(i) = a(i)
	    end do
		!$OMP END PARALLEL DO
		
	    if(ia(1) .eq. 1) then
			!$OMP PARALLEL DO PRIVATE(i)
			do i = 1, neqns+1
				ia(i) = ia(i)-1
			end do
			!$OMP END PARALLEL DO
			
			!$OMP PARALLEL DO PRIVATE(i)
			do i = 1, nnz
				ja(i) = ja(i)-1
			end do
			!$OMP END PARALLEL DO
		end if
	    
		call auction(ia, ja, tmp_a, neqns, nnz, match_perm, match_scaling, rBlockBndry, option, tot_neqns)
		
		if(rank .eq. 0) then
        	if(ia(1) .eq. 0) then
        		!$OMP PARALLEL DO PRIVATE(i)
        		do i = 1, neqns+1
        			ia(i) = ia(i)+1
        		end do
        		!$OMP END PARALLEL DO

        		!$OMP PARALLEL DO PRIVATE(i)
        		do i = 1, neqns
        			match_perm(i) = match_perm(i)+1
        		end do
        		!$OMP END PARALLEL DO
        		
        		!$OMP PARALLEL DO PRIVATE(i)
        		do i = 1, nnz
        			ja(i) = ja(i)+1
        		end do
        		!$OMP END PARALLEL DO
        	end if
		
		
            sumVal = 0.0d0
            do i = 1, neqns
               do j = ia(i), ia(i+1)-1
                 if(ja(j) .eq. match_perm(i)) then
                     sumVal = sumVal + abs(a(j))
                     exit
                  end if
               end do
    	     end do
             write(*,*) 'Weight of Auction Matching', sumVal
		end if
		
	end if

	deallocate(ia)
	deallocate(ja)
	deallocate(a)
	deallocate(tmp_a)

	if(rank .eq. 0) then
		deallocate(match_perm)
		deallocate(match_scaling)
		deallocate(rBlockBndry)
	end if

	call MPI_FINALIZE(error)

END PROGRAM main
