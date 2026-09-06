#ifndef CUDA_SPARSE_MATRIX_H
#define CUDA_SPARSE_MATRIX_H

#include <cusparse.h>
#include <thrust/device_vector.h>

#include "cuda_linalg/CudaUtil.h"

namespace cuda_linalg
{

/* CSR Sparse Matrix */
template<typename T>
class CudaSparseMatrix : public CudaObject
{
  public:

    CudaSparseMatrix()
      : num_col_     (0),
        num_row_     (0),
        mat_descr_   (nullptr),
        sparse_descr_(nullptr)
    {}

    int getNumNonzero() const { return static_cast<int>(d_values_.size()); }
    int getNumRow()     const { return num_row_; }
    int getNumCol()     const { return num_col_; }

    cusparseSpMatDescr_t getSparseDescriptor()       { return sparse_descr_; }
    cusparseSpMatDescr_t getSparseDescriptor() const { return sparse_descr_; }

    const CudaVector<T>& getValues() const { return d_values_; }
          CudaVector<T>& getValues()       { return d_values_; }

    const CudaVector<int>& getCRows () const { return d_row_comp_; }
          CudaVector<int>& getCRows ()       { return d_row_comp_; }

    const CudaVector<int>& getCols  () const { return d_col_; }
          CudaVector<int>& getCols  ()       { return d_col_; }

    T* getBufferPtr() { return d_buffer_mv_.data(); }

    int* getRowIndexCompPtr() { return d_row_comp_.data(); }
    int* getColIndexPtr()     { return d_col_.data();      }
    T*   getValuesPtr()       { return d_values_.data();   }
  
    const int* getRowIndexCompPtr() const { return d_row_comp_.data(); }
    const int* getColIndexPtr()     const { return d_col_.data();      }
    const T*   getValuesPtr()       const { return d_values_.data();   }
  
    void initialize(
      int num_col,
      const std::vector<int>& crows,
      const std::vector<int>& cols,
      const std::vector<T>&   nnzs)
    {
      cudaDataType_t data_type;
      if constexpr (std::is_same_v<T, float>)
        data_type = CUDA_R_32F;
      else if constexpr (std::is_same_v<T, double>)
        data_type = CUDA_R_64F;
      else
        assert(0);

      num_row_ = static_cast<int>(crows.size() - 1);
      num_col_ = num_col;
	    const int num_nonzero = static_cast<int>(nnzs.size());

      // Create Matrix Descriptor
      CHECK_CUSPARSE( cusparseCreateMatDescr(&mat_descr_) )

      // Reset and Copy Data
      d_row_comp_.clear();
      d_col_.clear();
      d_values_.clear();

      d_row_comp_.resize(num_row_ + 1);
      d_col_.resize(num_nonzero);
      d_values_.resize(num_nonzero);

    	std::vector<int> crow_index(crows.begin(), crows.end());
      std::vector<int> col_index(cols.begin(), cols.end());
      std::vector<T> values(nnzs.begin(), nnzs.end());

      // Copy Host to Device
      thrust::copy(crow_index.begin(),
                   crow_index.end(),
                   d_row_comp_.begin());
    
      thrust::copy(col_index.begin(),
                   col_index.end(),
                   d_col_.begin());
    
      thrust::copy(values.begin(),
                   values.end(),
                   d_values_.begin());

      CHECK_CUSPARSE( 
        cusparseCreateCsr(
          &sparse_descr_, /* spMatDescr : Sparse Matrix Descriptor */
          num_row_,  
          num_col_, 
          num_nonzero,
          (void*)(getRowIndexCompPtr()), /* int* must be casted to void* */
          (void*)(getColIndexPtr()),
          (void*)(getValuesPtr()),
          CUSPARSE_INDEX_32I,
          CUSPARSE_INDEX_32I,
          CUSPARSE_INDEX_BASE_ZERO, 
          data_type) )  /* cudaDataType */ 
    
      T alpha = 1.0;
      T beta  = 0.0;

      cusparseDnVecDescr_t vec_x_descr, vec_y_descr;

      CudaVector<T> d_vec_x(num_col_);
      CudaVector<T> d_vec_y(num_row_);
      T* d_vec_x_ptr = d_vec_x.data();
      T* d_vec_y_ptr = d_vec_y.data();

      CHECK_CUSPARSE( cusparseCreateDnVec(&vec_x_descr, num_col_, (void*)(d_vec_x_ptr), data_type) )
      CHECK_CUSPARSE( cusparseCreateDnVec(&vec_y_descr, num_row_, (void*)(d_vec_y_ptr), data_type) )

      size_t buffer_size_mv = 0;
      /* For NonTranspose(Matrix)-Vector Mult */
      CHECK_CUSPARSE(
        cusparseSpMV_bufferSize(
          getCuSparseHandle(),
          CUSPARSE_OPERATION_NON_TRANSPOSE,
          &alpha,
          sparse_descr_,
          vec_x_descr,
          &beta,
          vec_y_descr,
          data_type,
          CUSPARSE_SPMV_ALG_DEFAULT,
          &buffer_size_mv) )

      d_buffer_mv_.resize(buffer_size_mv);
      T* d_buffer_ptr = d_buffer_mv_.data();

      /* preprocess makes buffer 'active' and accelerates SpMV. This is optional. */
      CHECK_CUSPARSE(
        cusparseSpMV_preprocess(
          getCuSparseHandle(),
          CUSPARSE_OPERATION_NON_TRANSPOSE,
          &alpha,
          sparse_descr_,
          vec_x_descr,
          &beta,
          vec_y_descr,
          data_type,
          CUSPARSE_SPMV_ALG_DEFAULT,
          d_buffer_ptr) )

      CHECK_CUSPARSE( cusparseDestroyDnVec(vec_x_descr) )
      CHECK_CUSPARSE( cusparseDestroyDnVec(vec_y_descr) )
    }

  private:

    int num_row_;                         // NumRow
    int num_col_;                         // NumCol

    CudaVector<T>  d_buffer_mv_;          // Buffer for Matirx -Vector Mult

    cusparseMatDescr_t   mat_descr_;      // Matrix Descriptor
    cusparseSpMatDescr_t sparse_descr_;   // Sparse Matrix Descriptor
    CudaVector<int> d_row_comp_;          // Row    Index (Compressed)
    CudaVector<int> d_col_;               // Column Index
    CudaVector<T>   d_values_;            // NonZero Values
};

}

#endif
