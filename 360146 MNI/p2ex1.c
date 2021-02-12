#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void triang(int, double**, double*, int);

int main(void){
    int i, j, n, t;
    double **a, *b;
    
    printf("Tipus de sistemes a resoldre:\n\t0 -> triangular inferior\n\t1 -> triangular superior\n\t-1 -> triangular inferior amb 1 a la diagonal\nTipus: ");
    scanf("%d",&t);
    printf("Nombre d'equacions: ");
    scanf("%d", &n);
    
    a = (double**)malloc(n*sizeof(double*));
    if(a==NULL){
        printf("Memòria insuficient.\n");
        exit(1);
    }
    for (i=0; i<n; i++) {
        a[i] = (double*)malloc(n*sizeof(double));
        if(a[i]==NULL){
            printf("Memòria insuficient.\n");
            exit(2);
        }
    }
    b = (double*)malloc(n*sizeof(double));
    if(b==NULL){
        printf("Memòria insuficient.\n");
        exit(3);
    }
    
    printf("Matriu del sistema:\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("A[%d][%d] = ", i, j);
            scanf("%lf", &a[i][j]);
        }
    }
    
    printf("Termes independents:\n");
    for(i=0; i<n; i++){
        printf("b[%d] = ", i);
        scanf("%lf", &b[i]);
    }
    
    triang(n, a, b, t);
    
    printf("El vector solucio es:\n");
    for(i=0; i<n; i++){
        printf("x[%d] = %lf\n", i, b[i]);
    }
    
    for(i=0; i<n; i++){
        free(a[i]);
    }
    free(a);
    free(b);
    
    return 0;
}

void triang(int n, double **matA, double *b, int tipus){
    int i, j;
    double sum;
    switch (tipus) {
        case 1:
            for(i = 1; i <= n; i++){
                sum = 0;
                for(j = 1; j < i; j++){
                    sum = sum + matA[n-i][n-j]*b[n-j];
                }
                b[n-i] = (b[n-i]-sum)/(matA[n-i][n-i]);
            }
            break;
        case 0:
            for(i=n; i>=1; i--){
                sum = 0;
                for(j=n; j>i; j--){
                    sum = sum + matA[n-i][n-j]*b[n-j];
                }
                b[n-i] = (b[n-i]-sum)/matA[n-i][n-i];
            }
            break;
        case -1:
            for(i=n-1; i>=1; i--){
                sum = 0;
                for(j=n; j>i; j--){
                    sum = sum + matA[n-i][n-j]*b[n-j];
                }
                b[n-i] = (b[n-i]-sum);
            }
            break;
        default:
            printf("Tipus de matrius d'entrada incorrecte.");
            break;
    }
}
