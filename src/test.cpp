#include "iostream"
#include "Matrix.h"

using namespace std;
using namespace Matrix_ly;

int main()
{
    double a[8] = {1,2,3,4,5,6,7,8};
    double b[6] = {1,2,3,4,5,6};
    double c[25]= {4,2,5,6,7,
        2,3,4,2,1,
        7,8,2,3,5,
        6,4,2,4,6,
        9,0,8,7,3};
    //doublecomplex d[4] = {doublecomplex(2,0), doublecomplex(3,0),
    //  doublecomplex(8,9),doublecomplex(0,1)};

    //drvec aa(a,3);
    //dcvec bb(d,3);
    //cout << bb - aa << endl;

    double ti[9] = {5.29484e-17, -6.83136e-18, -3.87505e-18, -4.19491e-18, 5.81811e-17, -5.35767e-18, -2.80359e-18, -4.64279e-18, 6.57287e-17};
    drmat t(ti, 3, 3, false);
    cout << "ti = " << t << endl;

    //dcvec vlamda;
    //dcmat q;
    //dcmat mlamda;
    //dcmat qinv;
    //cout <<"eig success: " << (t.eig(vlamda,q)) << endl;
    //t.eig(vlamda,q);
    //cout << vlamda << endl << endl;
    //cout << q << endl << endl;
    //dcmat gg(d, 2, 2);
    //dcmat hh(d, 2, 2);
    //cout << gg << endl << endl;
    //cout << gg + hh << endl << endl;
    //dcvec vlamda;
    //dcmat q;
    //dcmat mlamda;
    //dcmat qinv;
    //cout <<"eig success: " << gg.eig(vlamda,q) << endl;
    //cout << vlamda << endl << endl;
    //cout << q << endl << endl;
    //qinv = q.inv();
    //mlamda.diag(vlamda); 
    //cout << q * mlamda * qinv << endl;
    //cout << "the root:" << endl;
    //dcmat qroot = gg.expm();
    //cout << qroot << endl << endl;
    //drmat mm(a, 2, 2);
    //dcmat mmroot;
    //mmroot = mm.expm();
    //cout << mmroot << endl << endl;
    drmat x1(c, 6, 4);
    drmat y1(a, 6, 1);

    cout << "a = " << x1 << endl << endl;
    cout << "b = " << y1 << endl << endl;
    cout << x1.backslash(y1) << endl << endl;

    //dcvec yy(d, 4);
    //cout << yy << endl << endl;

    //doublecomplex xx = 3.0;
    //cout << xx << endl << endl;


    //doublecomplex a(1.9, 2.3);
    //doublecomplex b(1.9, 2.3);
    //cout << a << endl;
    //cout << (a <= b) << endl;
    return 0;
}
