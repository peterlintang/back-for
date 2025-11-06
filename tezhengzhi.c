#include <stdio.h>
#include <math.h>
#include <complex.h>
/*在这里定义和修改尝试的各种参数。*/
#define MIN -10  /*尝试区间的下限*/
#define MAX 10  /*尝试区间的上限*/
#define STEP 0.001  /*尝试过程的步长*/
#define PRE 0.01  /*判断为零的精度*/
/*定义矩阵类型，其中C_MATRIX类型储存的元素均为复数。*/
typedef complex double C_MATRIX[20][20];
typedef double MATRIX[20][20];
/*定义复数类型的函数，用于计算行列式的值。*/
complex double determinant(C_MATRIX c_mat,int order);
complex double cofactor(C_MATRIX c_mat,int order,int r,int c);

int main(int argc,char *argv[])
{
    double real,imag;
    MATRIX mat;
    C_MATRIX c_mat,c_mat_alt;
    int order,i,j,sign = 0;
    
    printf("输入矩阵的阶数:");
    scanf("%d",&order);
    printf("输入矩阵:\n");
    for(i = 0;i < order;i ++)
        for(j = 0;j < order;j ++)
        {
            scanf("%lf",&mat[i][j]);
            /*将输入到mat中的元素赋值给存储复数类型元素的矩阵c_mat。*/
            c_mat[i][j] = (complex double)mat[i][j];
        }
    for(i = 0;i < order;i ++)
        for(j = 0;j < order;j ++)
        {
            printf("real: %lf image: %f\n", creal(c_mat[i][j]), cimag(c_mat[i][j]));
        }
    /*试根求复特征值。*/
    printf("特征值为:\n");
    for(real = MIN;real <= MAX;real += STEP)
        for(imag = MIN;imag <= MAX;imag += STEP)
        {
            /*以sign为0或1判断正的虚部之前是否需要输出加号。*/
            sign = 0;
            /*复制一份矩阵参与运算。*/
            for(i = 0;i < order;i ++)
                for(j = 0;j < order;j ++)
                    c_mat_alt[i][j] = c_mat[i][j];
            for(i = 0;i < order;i ++)
                c_mat_alt[i][i] -= (real + imag * I);
            /*如果行列式计算结果显示实部和虚部均小于给定精度（被认作为0），判定real + imag i是特征值。*/
            if(fabs(creal(determinant(c_mat_alt,order))) <= PRE && fabs(cimag(determinant(c_mat_alt,order))) <= PRE)
                {
                    /*调整输出格式，去掉不必要输出的部分。*/
                    if(fabs(real) > PRE)
                    {
                        printf("%.3f ",real);
                        sign = 1;
                    }
                    if(fabs(imag) > PRE)
                    {
                        if(imag > PRE)
                        {
                            if(sign == 1)
                                printf("+ ");
                            printf("%.3fi",imag);
                        }
                        if(imag < - PRE)
                            printf("%.3fi",imag);
                    }
                    if(fabs(real) < PRE && fabs(imag) < PRE)
                        printf("0.000");
                    printf("\n");
                }
        }
    
    return 0;
}
/*下面是求行列式的相关函数。*/
complex double determinant(C_MATRIX c_mat,int order)
{
    complex double result = 0;
    int i;
    
    if(order == 1)
        result = c_mat[0][0];
    else
        for(i = 0;i < order;i ++)
            result += pow(-1,i) * c_mat[i][0] * cofactor(c_mat,order,i,0);
    
    return result;
}

complex double cofactor(C_MATRIX c_mat,int order,int r,int c)
{
    complex double result = 0;
    C_MATRIX c_cofactor;
    int original_i,original_j,i,j;
    
    for(i = 0;i < order;i ++)
        for(j = 0;j < order;j ++)
        {
            original_i = i;
            original_j = j;
            if(i == r || j == c);
            else
            {
                if(i > r)
                    i --;
                if(j > c)
                    j --;
                c_cofactor[i][j] = c_mat[original_i][original_j];
                i = original_i;
                j = original_j;
            }
        }
    if(order >= 2)
        result = determinant(c_cofactor,order - 1);
    
    return result;
}

