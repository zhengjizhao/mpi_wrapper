!this program enables multiple instances of an application to run simultaneously.  
program mwrapper 

#ifdef INTEL
  use ifport  
#endif
  use mpi_wrapper
  implicit none
  include 'mpif.h'
  integer :: global_comm = MPI_COMM_WORLD, global_size, global_rank, njobs, color, key, size, rank, ierr
  integer, parameter :: strlen = 200  
  integer ::  ipartition = 0  
  character(len=strlen) :: mydir, stdoutfile, stderrfile 
#ifdef  SERIAL
  character(len=strlen) :: appname 
#endif
  double precision ::  tstart, tend 

! initial setup
  call MPI_Init(ierr)
  call MPI_Comm_rank(global_comm, global_rank, ierr)
  call MPI_Comm_size(global_comm, global_size, ierr)

! process inputs, stop unpon error
#ifdef SERIAL 
  call process_input(global_comm, global_size, global_rank, njobs, appname, mydir, *500)
#else
  call process_input(global_comm, global_size, global_rank, njobs, mydir, *500)
#endif

! split processors into "njobs" number of groups, "color" is the gorup index
  call groups(global_rank, global_size, njobs, color, key, ipartition)
  call MPI_Comm_split(global_comm, color, key, aMPI_COMM_WORLD,ierr)
  call MPI_Comm_rank(aMPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(aMPI_COMM_WORLD, size, ierr)

  call print_msg(global_rank, 'njobs, global_size, global_size / njobs', 400)
#ifdef SERIAL  
  call print_msg(global_rank, 'application is ' // trim(appname), 499)
#endif

  call print_msg(rank, 'jobs starting', 401)

! cd to run directory 
  ierr = chdir(mydir)

! redirect the standard output (unit=6) and error (unit-0) to files for each job 
  write(stdoutfile, "('stdoutfile.', i4.4)") color
  write(stderrfile, "('stderrfile.', i4.4)") color
  if (rank == 0) open(unit = 6,file = stdoutfile, action = 'write')
  if (rank == 0) open(unit = 0,file = stderrfile, action = 'write')
  
! run njobs number of jobs simultantously, one per group
  tstart = MPI_Wtime()
#ifdef SERIAL 
  call execute_command_line(appname, wait = .true., exitstat = ierror_app)
#else
!  call MPI_Comm_set_errhandler(aMPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr) 
  call vamp()  !can be any mpi application (need to convert to a subroutine first)  
!  call jacobi_mpiomp() !this is an example for jacobi code
#endif
  tend = MPI_Wtime()

  if (rank == 0) close(unit = 6)
  if (rank == 0) close(unit = 0)
#ifdef GNU
  if (rank == 0) open(unit = 6, file = "/dev/stdout",status='old')  !reconnect to stdout device
#endif
  call print_msg(rank, 'jobs ended', 402)

!  call MPI_Barrier(global_comm,ierr)

! finish up
500  call MPI_Finalize(ierr)

  contains 

! only rank 0 prints
  subroutine print_msg(myrank, message, label)
  integer, intent(in) :: myrank
  character(*), intent(in), optional :: message
  integer, intent(in) :: label
  character(len=50) :: job_status

400   format('running ', I4, ' jobs with ', I10, ' processors, ', I4, ' processors each.')
401   format('the ', I4, '-th job with ', I4, ' processors running in directory, ', A, ' ...')  
402   format('the ', I4, '-th job running in ', A,  ' is ', A, '. elapsed time (sec): ', F10.2)  

  if (myrank == 0) then
      select case (label)
          case(400)
              write(*, 400) njobs, global_size, global_size/njobs 
          case(401)
              write(*, 401) color, size, trim(mydir) 
          case(402) 
              job_status = 'completed'
              if (ierror_app /= 0) job_status = 'failed'
              write(*, 402) color, trim(mydir), trim(job_status), tend - tstart 
          case default
              write(*, *) message
      end select
  endif

  return
  end subroutine

end program

#ifdef SERIAL 
  subroutine process_input(global_comm, global_size, global_rank, njobs, appname, mydir, *)
#else
  subroutine process_input(global_comm, global_size, global_rank, njobs, mydir, *)
#endif
  include 'mpif.h'
  integer, intent(in) :: global_comm, global_size, global_rank
  integer, intent(out) :: njobs
  integer, parameter :: strlen = 200
#ifdef SERIAL 
  character(len = strlen), intent(out) :: appname 
#endif
  character(len = strlen), intent(out) :: mydir 
  character(len = strlen), allocatable, dimension(:) :: rundirs
  character(len = strlen) :: joblistfile
  integer :: ierr, ierror, nlines, i
  logical :: iexist

! only the global rank 0 read in the joblist file
  ierror = 0 
  if (global_rank == 0) then

!     defaults  
      joblistfile = 'joblist.in'
#ifdef SERIAL
      appname = 'myapp.x'
#endif
!     get optional command line arguments
      nargs = command_argument_count()
      if (nargs >= 1)  call get_command_argument(1,joblistfile)
#ifdef SERIAL
      if (nargs >= 2)  call get_command_argument(2,appname)
#endif

      inquire(file = joblistfile, exist = iexist)
      if (.not. iexist) then
           print *, 'the joblistfile ', trim(joblistfile), ' does not exist, stop'
           !call MPI_Abort(global_comm, 1, ierr) 
           ierror = 1
           goto 200
      endif

#ifdef SERIAL
      call find_app(appname, ierr)
      if (ierr /= 0) then
           print *, 'the executable ', trim(appname), ' does not exist, stop'
           !call MPI_Abort(global_comm, 1, ierr) 
           ierror = 1
           goto 200
      endif
#endif

      open (10,file = trim(joblistfile), status = 'old', action = 'read',iostat = ierr)
      rewind(10)
      read(10,*) njobs
      allocate(rundirs(0:njobs-1))
      nlines = 0
      do i = 0, njobs - 1
          read(10, '(A)', END = 100) rundirs(i)
          nlines = nlines + 1
      enddo

100   continue
      close(10)

!     sanity checking:
      if (nlines /= njobs) then
           print *, 'the number of run diretories does not match the number of jobs provided &
           &            on the first line of joblistfile, stop'
           print *, 'nlines = ', nlines, ' njobs = ', njobs
           !call MPI_Abort(global_comm, 1, ierr) 
           ierror = 1
           goto 200
      endif

      ncores = global_size / njobs
#ifdef SERIAL
      if (ncores /= 1) then 
           print *, 'more than one tasks are assigned to a serial execution, stop' 
           print *, 'ncores (tasks per job) =', ncores  
           ierror = 1
           goto 200
      endif
#endif
      if ( ncores * njobs /= global_size ) then
           print *, 'the total number of processors is not divisble by the number of jobs, stop ' 
           print *, 'global_size, njobs = ', global_size, njobs  
           !call MPI_Abort(global_comm, 1, ierr) 
           ierror = 1
           goto 200
      endif

      do i = 0, njobs - 1
#ifdef INTEL
          inquire(directory = rundirs(i), exist = iexist) 
#else  
          inquire(file = rundirs(i), exist = iexist)
#endif
          if (iexist .eqv. .false.) then
               print *, 'the run directory ', trim(rundirs(i)), ' does not exit, stop'
               !call MPI_Abort(global_comm, 1, ierr) 
               ierror = 1
               goto 200
          endif
      enddo
  endif

200 call MPI_Bcast(ierror, 1, MPI_integer, 0, global_comm, ierr)

  if (ierror /= 0 ) return 1  !stop upon error 

! broadcast jobinfo for the rest of the ranks 
  call MPI_Bcast(njobs, 1, MPI_integer, 0, global_comm, ierr)
#ifdef SERIAL  
  call MPI_Bcast(appname, strlen, MPI_character, 0, global_comm, ierr)
#endif
  if (global_rank == 0) then
      mydir = rundirs(0)
      do i = 1, global_size - 1
         ii = i / ncores  !equivalent to color 
         call MPI_Send(rundirs(ii), strlen, MPI_Character, i, i, global_comm,ierr)
      enddo
  else
      call MPI_Recv(mydir, strlen, MPI_Character, 0, global_rank, global_comm, status, ierr)
  endif

  if (global_rank == 0) then
      deallocate(rundirs)
  endif

  return

end subroutine

#ifdef SERIAL
! get the absolute path to the application
  subroutine find_app(app, ierr)
      character(*), intent(inout) :: app
      integer, intent(out) :: ierr

      call execute_command_line('which ' // app // '> /tmp/fullpath.file 2>/dev/null', wait = .true., exitstat = ierr)
      if (ierr == 0) then
          open(unit = 12, file = '/tmp/fullpath.file', status = 'old', action = 'read', iostat = ierr)
          read(12, "(A)") app  ! "(A)" works but not "*"
          close(12)
          call execute_command_line('rm -fr /tmp/fullpath.file', wait=.false.)
      endif

      return
  end subroutine
#endif

! different ways of splitting processors into sub groups can be provided here 
subroutine groups(global_rank, global_size, njobs, color, key, ischeme)
    integer, intent(in) :: global_rank, global_size, njobs, ischeme
    integer, intent(out) :: color, key 
    integer :: ncores

    if (ischeme == 0) then
        ncores = global_size / njobs
        color = global_rank / ncores 
        key = global_rank - color * ncores
    else 
        color = mod(global_rank, njobs)
        key = int(global_rank / njobs)
    endif

    return
end subroutine groups 

