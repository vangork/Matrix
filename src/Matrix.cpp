#include "Matrix.h"

namespace Matrix_ly
{
    void dcmat2dr(const dcmat &__m, drmat &__n)
    {
        int _row = __m.row();
        int _col = __m.col();
        //assert(_row == __n.row() && "Matrix rows must be equal");
        //assert(_col == __n.col() && "Matrix columns must be equal");
        __n.set(_row, _col);
        for(int _i = 0; _i < _row * _col; ++_i)
        {
            *(__n.data() + _i) = (__m.data() + _i)->real();
        }
    }

    void dcvec2dr(const dcvec &__m, drvec &__n)
    {
        int _len = __m.len();
        //assert(_len == __n.len() && "Vectors len must be equal");
        __n.set(_len);
        for(int _i = 0; _i < _len; ++_i)
        {
            *(__n.data() + _i) = (__m.data() + _i)->real();
        }
    }

}
