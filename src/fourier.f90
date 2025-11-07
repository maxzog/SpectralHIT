module spectral
   use, intrinsic :: ISO_C_BINDING
   ! implicit none
   include 'fftw3.f03'

   ! Utility
   real(8), parameter :: PI=3.141592653589793
   double complex, parameter :: i_unit=cmplx(0., 1.)

   contains

      ! Return centered wavenumber vector from 0 to N/2-1, -N/2 to -1 assuming domain in (0, 2pi)
      function get_kv(rank) result(kv)
         implicit none
         integer, dimension(rank) :: kv
         integer :: rank, i
         do i=1, rank
            if (i <= rank/2) then
               kv(i) = i - 1
            else
               kv(i) = i - rank - 1
            end if
         end do
      end function get_kv

      ! Compute the autocorrelation of arr in Fourier space
      ! assumes the original function is real
      function compute_autocorr(arr, rank) result(out)
         implicit none
         double complex, dimension(rank) :: arr, out
         integer :: rank, i
         do i=1,rank
            out(i) = arr(i)*conjg(arr(i))
         end do
      end function compute_autocorr

      ! Compute the energy spectrum of arr
      function compute_spectrum(arr, rank) result(out)
         implicit none
         double complex, dimension(rank) :: arr, out
         integer :: rank, i
         do i=1,rank
            out(i) = arr(i)*conjg(arr(i))
         end do
      end function compute_spectrum

      ! Compute the radial energy spectrum of 3D arr
      subroutine compute_radial_spectrum(arr, mesh, out, axis)
         use grid_class, only: grid
         implicit none
         class(grid), intent(in) :: mesh
         double complex, dimension(:,:,:,:), intent(in) :: arr
         real(8), dimension(:), intent(inout) :: out, axis
         integer, dimension(:), allocatable :: counts
         real(8) :: kv(3),dk,kmag
         integer :: i,j,k,nbins,ind

         out=0.0d0
         nbins=size(out)

         ! Allocate counter 
         allocate(counts(1:nbins)); counts=0

         ! Compute bin properties
         kv=[real(mesh%nx/2+1, 8), real(mesh%ny/2+1,8), real(mesh%nz/2+1,8)]
         kmag=norm2(kv)
         dk=kmag/nbins

         ! Compute spectrum
         do k=1,mesh%Nz
            do j=1,mesh%Ny
               do i=1,mesh%Nx
                  kv = [mesh%k1v(i),mesh%k2v(j),mesh%k3v(k)]
                  kmag = norm2(kv)
                  if (kmag.lt.10.0d0*epsilon(1.0d0)) cycle
                  ind = ceiling(kmag/dk)
                  out(ind) = out(ind) + REAL(arr(1,i,j,k)*conjg(arr(1,i,j,k))) & 
                           &          + REAL(arr(2,i,j,k)*conjg(arr(2,i,j,k))) & 
                           &          + REAL(arr(3,i,j,k)*conjg(arr(3,i,j,k)))
                  counts(ind) = counts(ind) + 1
               end do
            end do
         end do

         ! Create axis (for plotting) 
         do i=1,nbins
            axis(i)=dk*(i-1)
         end do

         deallocate(counts)
      end subroutine compute_radial_spectrum

      !> ONE-DIMENSIONAL ROUTINES
      ! Subroutine that computes the one-dimensional Fourier transform
      subroutine FFT_1D(in, out, rank0)
         implicit none
         integer(C_INT), intent(in) :: rank0                        ! Number of wavenumbers
         double complex, dimension(rank0), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_1d(rank0, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, in, out)
         ! Normalize
         out = out/rank0
         ! Cleanup
         call fftw_destroy_plan(plan)
      end subroutine FFT_1D

      ! Subroutine that computes the one-dimensional Fourier transform
      subroutine iFFT_1D(in, out, rank0)
         implicit none
         integer(C_INT), intent(in) :: rank0                        ! Number of wavenumbers
         double complex, dimension(rank0), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_1d(rank0, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, out, in)
         ! Cleanup
         call fftw_destroy_plan(plan)
      end subroutine iFFT_1D

      subroutine ddx_1D(in, wv, rank0)
         implicit none
         integer, intent(in) :: rank0
         double complex, dimension(rank0), intent(inout) :: in
         integer, dimension(rank0), intent(in) :: wv
         integer :: i

         do i=1,rank0
            in(i) = wv(i)*i_unit*in(i)
         end do
      end subroutine ddx_1D

      subroutine initialize_1D(in, wv, rank0)
         implicit none
         double complex, dimension(rank0), intent(inout) :: in
         integer, dimension(rank0), intent(in) :: wv
         integer, intent(in) :: rank0
         integer :: i, k

         ! Set oddball element
         in(rank0/2+1)=cmplx(1.,0.)
         ! Set zeroth mode
         in(1)=cmplx(1.,0.)

         ! Fill in the rest from the outside in
         do i=2,rank0/2
            k = wv(i)
            in(rank0-(i-2))=cmplx(real(i,8), k, kind=8)
            in(i)=conjg(in(rank0-(i-2)))
         end do
      end subroutine initialize_1D

      !> TWO-DIMENSIONAL ROUTINES
      ! Subroutine that computes the two-dimensional Fourier transform
      subroutine FFT_2D(in, out, rank0, rank1)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1                        ! Number of wavenumbers
         double complex, dimension(rank0, rank1), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_2d(rank1, rank0, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, in, out)
         ! Normalize
         out = out/(rank0*rank1)
      end subroutine FFT_2D

      ! Subroutine that computes the two-dimensional Fourier transform
      subroutine iFFT_2D(in, out, rank0, rank1)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1                        ! Number of wavenumbers
         double complex, dimension(rank0, rank1), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_2d(rank1, rank0, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, out, in)
      end subroutine iFFT_2D

      subroutine initialize_2D(in, w1v, w2v, rank0, rank1)
         implicit none
         double complex, dimension(rank0, rank1), intent(inout) :: in  ! 2D complex array
         integer, dimension(rank0), intent(in) :: w1v                  ! Wavenumbers in x direction
         integer, dimension(rank1), intent(in) :: w2v                  ! Wavenumbers in y direction
         integer, intent(in) :: rank0, rank1                           ! Dimensions of the array
         integer :: i, j                                               ! Loop indices
         real(8) :: kmag
         in=cmplx(0.,0.)

         ! Set zeroth mode and oddball values
         in(1,1)=cmplx(1.,0.)
         in(rank0/2+1,1)=cmplx(1.,0.)
         in(1,rank1/2+1)=cmplx(1.,0.)

         ! Fill in the rest from the outside in
         do i=2,rank0/2
            do j=2,rank1/2
               kmag = sqrt(real(w1v(i)**2 + w2v(j)**2, 8))
               in(rank0-(i-2),rank1-(j-2))=cmplx(real(i,8), kmag, kind=8)
               in(i,j)=conjg(in(rank0-(i-2), rank1-(j-2)))
            end do
         end do
      end subroutine initialize_2D

      ! Subroutine that computes the three-dimensional Fourier transform
      subroutine FFT_3D(in, out, rank0, rank1, rank2)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1, rank2                        ! Number of wavenumbers
         double complex, dimension(rank0, rank1, rank2), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_3d(rank2, rank1, rank0, in, out, FFTW_FORWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, in, out)
         ! Normalize
         out = out/(rank0*rank1*rank2)
         ! Cleanup
         call fftw_destroy_plan(plan)
      end subroutine FFT_3D

      ! Subroutine that computes the three-dimensional Fourier transform
      subroutine iFFT_3D(in, out, rank0, rank1, rank2)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1, rank2                        ! Number of wavenumbers
         double complex, dimension(rank0, rank1, rank2), intent(inout) :: in, out ! Function and fourier coefficients
         type(C_PTR) :: plan
         
         ! Plan 
         plan=fftw_plan_dft_3d(rank2, rank1, rank0, out, in, FFTW_BACKWARD, FFTW_ESTIMATE)
         ! Execute
         call fftw_execute_dft(plan, out, in)
         ! Cleanup
         call fftw_destroy_plan(plan)
      end subroutine iFFT_3D

      subroutine initialize_3D(in, w1v, w2v, w3v, rank0, rank1, rank2)
         implicit none
         double complex, dimension(rank0, rank1, rank2), intent(inout) :: in  ! 2D complex array
         integer, dimension(rank0), intent(in) :: w1v                         ! Wavenumbers in x direction
         integer, dimension(rank1), intent(in) :: w2v                         ! Wavenumbers in y direction
         integer, dimension(rank2), intent(in) :: w3v                         ! Wavenumbers in y direction
         integer, intent(in) :: rank0, rank1, rank2                           ! Dimensions of the array
         integer :: i, j, k                                                   ! Loop indices
         real(8) :: kmag                                                      ! Wavevector 2-norm
         in=cmplx(0.,0.)

         ! Set zeroth mode and oddball values
         in(1,1,1)=cmplx(1.,0.)
         in(rank0/2+1,1,1)=cmplx(1.,0.)
         in(1,rank1/2+1,1)=cmplx(1.,0.)
         in(1,1,rank2/2+1)=cmplx(1.,0.)

         ! Fill in the rest from the outside in
         do k=2,rank2/2
            do j=2,rank1/2
               do i=2,rank0/2
                  kmag = sqrt(real(w1v(i)**2 + w2v(j)**2 + w3v(k)**2,8))
                  in(rank0-(i-2),rank1-(j-2),rank2-(k-2))=cmplx(real(i,8), kmag, kind=8)
                  in(i,j,k)=conjg(in(rank0-(i-2), rank1-(j-2), rank2-(k-2)))
               end do
            end do
         end do
      end subroutine initialize_3D

      ! Routine to write complex array (1D of size N) to text file
      subroutine write_complex_array_1D(arr, rank0, filename)
         implicit none
         integer(C_INT), intent(in) :: rank0            
         double complex, dimension(rank0), intent(in) :: arr 
         character(len=*), intent(in) :: filename  
         integer :: line
         real(8) :: real_part, imag_part
         integer :: unit
      
         unit = 1  
         open(unit, file=filename, status='replace', action='write', form='formatted')
      
         do line = 1, rank0
             real_part = real(arr(line))  
             imag_part = aimag(arr(line)) 
             write(unit, '(F16.8, 1X, F16.8)') real_part, imag_part
         end do
      
         close(unit)
      end subroutine write_complex_array_1D

      ! Routine to write complex array (2D of size rank0 x rank1) to text file
      subroutine write_complex_array_2D(arr, rank0, rank1, filename)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1                   ! Dimensions of the array
         double complex, dimension(rank0, rank1), intent(in) :: arr   ! 2D complex array 
         character(len=*), intent(in) :: filename                     ! Output filename
         integer :: i, j                                              ! Loop variables for 2D iteration
         real(8) :: real_part, imag_part                              ! Variables for real and imaginary parts
         integer :: unit                                              ! Unit number for the file
      
         unit = 1  
         open(unit, file=filename, status='replace', action='write', form='formatted')
      
         ! Iterate over both dimensions of the 2D array and write real and imaginary parts to file
         do j = 1, rank1
            do i = 1, rank0
               real_part = real(arr(i,j))  
               imag_part = aimag(arr(i,j)) 
               write(unit, '(F16.8, 1X, F16.8)') real_part, imag_part
            end do
         end do
      
         close(unit)
      end subroutine write_complex_array_2D

      ! Routine to write complex array (3D of size rank0 x rank1 x rank2) to text file
      subroutine write_complex_array_3D(arr, rank0, rank1, rank2, filename)
         implicit none
         integer(C_INT), intent(in) :: rank0, rank1, rank2            ! Dimensions of the array
         double complex, dimension(rank0, rank1, rank2), intent(in) :: arr   ! 2D complex array 
         character(len=*), intent(in) :: filename                     ! Output filename
         integer :: i, j, k                                           ! Loop variables for 2D iteration
         real(8) :: real_part, imag_part                              ! Variables for real and imaginary parts
         integer :: unit                                              ! Unit number for the file
      
         unit = 1  
         open(unit, file=filename, status='replace', action='write', form='formatted')
      
         ! Iterate over both dimensions of the 2D array and write real and imaginary parts to file
         do k = 1, rank2
            do j = 1, rank1
               do i = 1, rank0
                  real_part = real(arr(i,j,k))  
                  imag_part = aimag(arr(i,j,k)) 
                  write(unit, '(F16.8, 1X, F16.8)') real_part, imag_part
               end do
            end do
         end do
         close(unit)
      end subroutine write_complex_array_3D

      ! Routine to write real array (1D of size N) to text file
      subroutine write_real_array(arr, rank0, filename)
         implicit none
         integer(C_INT), intent(in) :: rank0            
         double complex, dimension(rank0), intent(in) :: arr 
         character(len=*), intent(in) :: filename  
         real(8) :: val
         integer :: line
         integer :: unit
      
         unit = 1  
         open(unit, file=filename, status='replace', action='write', form='formatted')
      
         do line = 1, rank0
           val = real(arr(line))
           write(unit, '(F16.8)') val 
         end do
      
         close(unit)
      end subroutine write_real_array

      ! Routine to write real array (1D of size N) to text file
      subroutine write_double_array(arr, rank0, filename)
         implicit none
         integer(C_INT), intent(in) :: rank0            
         real(8), dimension(rank0), intent(in) :: arr 
         character(len=*), intent(in) :: filename  
         real(8) :: val
         integer :: line
         integer :: unit
      
         unit = 1  
         open(unit, file=filename, status='replace', action='write', form='formatted')
      
         do line = 1, rank0
           val = arr(line)
           write(unit, '(F16.8)') val 
         end do
      
         close(unit)
      end subroutine write_double_array
end module spectral
