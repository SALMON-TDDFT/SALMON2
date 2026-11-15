module cublas_wrapper
  use iso_c_binding
  use cublas
  implicit none

  type(cublasHandle), save :: cublas_handle
  logical, save :: handle_initialized = .false.

contains

  subroutine init_handle_if_needed
    use openacc
    implicit none

    integer(4) :: istat

    if (.not. handle_initialized) then
      istat = cublasCreate(cublas_handle)
      istat = cublasSetStream(cublas_handle, acc_get_cuda_stream(acc_async_sync))

      handle_initialized = .true.
    end if
  end subroutine init_handle_if_needed

  integer(4) function cublasZgemmStridedBatch_wrapper(transa, transb, m, n, k, alpha, Aarray, lda, strideA, Barray, ldb, strideB, beta, Carray, ldc, strideC, batchCount)
    implicit none

    integer, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k
    complex(8), intent(in) :: alpha, beta
    complex(8), device :: Aarray(*), Barray(*)
    complex(8), device :: Carray(*)
    integer, intent(in) :: lda, ldb, ldc
    integer, intent(in) :: strideA, strideB, strideC
    integer, intent(in) :: batchCount

    call init_handle_if_needed()

    cublasZgemmStridedBatch_wrapper = cublasZgemmStridedBatched(cublas_handle, transa, transb, m, n, k, alpha, Aarray, lda, strideA, &
      Barray, ldb, strideB, beta, Carray, ldc, strideC, batchCount)
  end function cublasZgemmStridedBatch_wrapper
end module
