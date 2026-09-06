#ifndef CUDA_MATRIX_H
#define CUDA_MATRIX_H

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cusparse.h>

#include "cuda_linalg/CudaUtil.h"
#include "cuda_linalg/CudaVector.h"

namespace cuda_linalg
{

/* Row-Wise Flatten */
template<typename T>
class CudaFlattenMatrix : public CudaObject
{
  public:

    CudaFlattenMatrix() : num_row_(0), num_col_(0) {}

    CudaFlattenMatrix(int num_row, int num_col)
      : num_row_(num_row),
        num_col_(num_col)
    {
      initialize(num_row, num_col);
    }

    int getNumEntry() const { return num_row_ * num_col_; }
    int getNumRow()   const { return num_row_; }
    int getNumCol()   const { return num_col_; }

          CudaVector<T>& getFlattenVector()       { return d_vector_; }
    const CudaVector<T>& getFlattenVector() const { return d_vector_; }

    cusparseDnVecDescr_t getVectorDescriptor()       { return d_vector_.getDescriptor(); }
    cusparseDnVecDescr_t getVectorDescriptor() const { return d_vector_.getDescriptor(); }

    void fillZero() { d_vector_.fillZero(); }

    void initialize(int num_row, int num_col)
    {
      d_vector_.clear();
      num_row_ = num_row;
      num_col_ = num_col;
      d_vector_.resize(num_row * num_col);
      d_vector_.fillZero();
    }

    void swap(CudaFlattenMatrix<T>& mat)
    {
      int mat_num_row = mat.getNumRow();
      int mat_num_col = mat.getNumCol();

      assert(mat_num_row == num_row_);
      assert(mat_num_col == num_col_);

      auto& d_mat_data = mat.getFlattenVector();

      d_vector_.swap(d_mat_data);
    }

  private:

    int num_row_;
    int num_col_;
    CudaVector<T> d_vector_;
};

}

#endif
