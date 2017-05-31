// File: array3d.h
// Author:Joe Barker
// Last-modified: 20 May 2014 12:22:20

#ifndef __ARRAY3D_H__
#define __ARRAY3D_H__

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#include <vector>
#include <cassert>
#include <cstddef>
#include <fftw3.h>

template <typename _Tp>
class Array3D
{
  public:
    typedef unsigned int size_type;
    typedef _Tp value_type;
    typedef _Tp* iterator;
    typedef const _Tp* const_iterator;
    typedef ptrdiff_t difference_type;

    Array3D() : dim0(0), dim1(0), dim2(0), data(0) {}

    Array3D(size_type d0, size_type d1, size_type d2)
      : dim0(d0), dim1(d1), dim2(d2), data(d0*d1*d2) {}

    ~Array3D(){data.clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      dim2 = 0;
      data.clear();
    }

    ///IFill means initialize or fill
    inline void IFill(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k]=a;
                }
            }
        }
    }
    inline _Tp* RESTRICT ptr() {
      return &(data[0]);
    }
    inline int getarrayelement(unsigned int i,unsigned int j,unsigned int k)
    {
        return((i*dim1+j)*dim2+k);
    }
    inline void resize(size_type d0, size_type d1, size_type d2) {
      dim0 = d0; dim1 = d1; dim2 = d2;
      data.resize(d0*d1*d2);
    }

    inline _Tp& RESTRICT operator()(const size_type i, const size_type j,
        const size_type k) {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

    inline const _Tp& operator()(const size_type i, const size_type j,
        const size_type k) const {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

//    inline size_type size(const size_type i) const { return dim[i]; }

  private:
    size_type dim0;
    size_type dim1;
    size_type dim2;
    std::vector<_Tp> data;
};


// template specialization for fftw
template <>
class Array3D<fftw_complex>
{
  public:
    typedef unsigned int size_type;
    typedef fftw_complex value_type;
    typedef fftw_complex* iterator;
    typedef const fftw_complex* const_iterator;
    typedef ptrdiff_t difference_type;

    Array3D() : dim0(0), dim1(0), dim2(0), data(NULL) {}

    Array3D(size_type d0, size_type d1, size_type d2)
      {resize(d0,d1,d2);}

    Array3D(size_type d0, size_type d1, size_type d2, const fftw_complex &init)
		{
			resize(d0,d1,d2);
			for(size_type i=0; i<(d0*d1*d2); ++i){
				data[i][0] = init[0];
				data[i][1] = init[1];
			}
		}

    ~Array3D(){clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      dim2 = 0;
	  fftw_free(data);
      data = NULL;
    }

    inline void resize(size_type d0, size_type d1, size_type d2) {
      dim0 = d0; dim1 = d1; dim2 = d2;
      if(data != NULL){
          fftw_free(data);
      }
	  data = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*d0*d1*d2);
    }
    ///IFill means initialize or fill
    inline void IFill(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][0]=a;
                    data[(i*dim1+j)*dim2+k][1]=a;
                }
            }
        }
    }
    inline int getarrayelement(unsigned int i,unsigned int j,unsigned int k)
    {
        return((i*dim1+j)*dim2+k);
    }

    ///IFill means initialize or fill the complex part only
    inline void IFillComplex(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][1]=a;
                }
            }
        }
    }
    ///IFill means initialize or fill the real part only
    inline void IFillReal(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][0]=a;
                }
            }
        }
    }



    inline fftw_complex& RESTRICT operator()(const size_type i, const size_type j,
        const size_type k) {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

    inline const fftw_complex& operator()(const size_type i, const size_type j,
        const size_type k) const {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

	inline fftw_complex* __restrict__ ptr() {
      return data;
	}

//    inline size_type size(const size_type i) const { return dim[i]; }

  private:
    size_type dim0;
    size_type dim1;
    size_type dim2;
    fftw_complex* data;
};

// template specialization for fftw float
template <>
class Array3D<fftwf_complex>
{
  public:
    typedef unsigned int size_type;
    typedef fftwf_complex value_type;
    typedef fftwf_complex* iterator;
    typedef const fftwf_complex* const_iterator;
    typedef ptrdiff_t difference_type;

    Array3D() : dim0(0), dim1(0), dim2(0), data(NULL) {}

    Array3D(size_type d0, size_type d1, size_type d2)
      {resize(d0,d1,d2);}

    Array3D(size_type d0, size_type d1, size_type d2, const fftwf_complex &init)
		{
			resize(d0,d1,d2);
			for(size_type i=0; i<(d0*d1*d2); ++i){
				data[i][0] = init[0];
				data[i][1] = init[1];
			}
		}

    ~Array3D(){clear();}

    inline void clear() {
      dim0 = 0;
      dim1 = 0;
      dim2 = 0;
	  fftwf_free(data);
      data = NULL;
    }

    inline void resize(size_type d0, size_type d1, size_type d2) {
      dim0 = d0; dim1 = d1; dim2 = d2;
      if(data != NULL){
          fftwf_free(data);
      }
	  data = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*d0*d1*d2);
    }
    ///IFill means initialize or fill
    inline void IFill(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][0]=a;
                    data[(i*dim1+j)*dim2+k][1]=a;
                }
            }
        }
    }
    inline int getarrayelement(unsigned int i,unsigned int j,unsigned int k)
    {
        return((i*dim1+j)*dim2+k);
    }

    ///IFill means initialize or fill the complex part only
    inline void IFillComplex(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][1]=a;
                }
            }
        }
    }
    ///IFill means initialize or fill the real part only
    inline void IFillReal(size_type a)
    {
        for(unsigned int i = 0 ; i < dim0 ; i++)
        {
            for(unsigned int j = 0 ; j < dim1 ; j++)
            {
                for(unsigned int k = 0 ; k < dim2 ; k++)
                {
                    data[(i*dim1+j)*dim2+k][0]=a;
                }
            }
        }
    }



    inline fftwf_complex& RESTRICT operator()(const size_type i, const size_type j,
        const size_type k) {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

    inline const fftwf_complex& operator()(const size_type i, const size_type j,
        const size_type k) const {
      assert( i < dim0 );
      assert( j < dim1 );
      assert( k < dim2 );
      return data[(i*dim1+j)*dim2+k];
    }

	inline fftwf_complex* __restrict__ ptr() {
      return data;
	}

//    inline size_type size(const size_type i) const { return dim[i]; }

  private:
    size_type dim0;
    size_type dim1;
    size_type dim2;
    fftwf_complex* data;
};

#endif // __ARRAY3D_H__
