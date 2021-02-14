#include<stdio.h>
#include<stdlib.h>
#include<math.h>

void triang(int, double**, double*, int);
int lu(double**, int, int*, double);

int main(void){
    int i, j;
    int n, m, nP, supInd=0;
    int *perm;
    FILE *fptr, *sort;
    double supMatSist, supMatSistAux=0, aux;
    double tol=1e-8;
    double *supMatSol, *termInd, *auxV;
    double **a;
    
    fptr = fopen("matriuSistema.txt", "r");
    if(fptr==NULL){
        printf("Error d'obertura de fitxer.");
        exit(1);
    }
    
    fscanf(fptr, "%d", &n);
    
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
    
    perm = (int*)malloc(n*sizeof(int));
    if(perm==NULL){
        printf("Memoria insuficient.\n");
        exit(1);
    }
    for(i=0; i<n; i++){
        perm[i] = i;
    }
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            fscanf(fptr, "%lf", &a[i][j]);
        }
    }
    
    printf("Matriu d'ordre %d d'entrada:\n", n);
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            printf("%14.9e  ", a[i][j]);
        }
        printf("\n");
    }
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            supMatSistAux += fabs(a[i][j]);
        }
        if(supMatSistAux>supMatSist){
            supMatSist = supMatSistAux;
        }
        supMatSistAux=0;
    }
    
    nP = lu(a, n, perm, tol);
    
    printf("Matriu L:\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i==j){
                printf("%12.9e  ", 1.0);
            }else if(i<j){
                printf("%12.9e  ", 0.0);
            }else{
                printf("%12.9e  ", a[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");
    
    printf("Matriu U:\n");
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(i>j){
                printf("%12.9e  ", 0.0);
            }else{
                printf("%12.9e  ", a[i][j]);
            }
        }
        printf("\n");
    }
    printf("\n");
    
    if(nP==-1){
        printf("No s'ha pogut fer la deescomposició.\n");
        exit(3);
    }
    
    sort = fopen("sistemesSolucionats.txt", "w");
    if(sort==NULL){
        printf("Error d'obertura de fitxer.");
        exit(1);
    }
    
    termInd = (double*)malloc(n*sizeof(double));
    if(termInd==NULL){
        printf("Memoria insuficient.\n");
        exit(1);
    }
    
    fscanf(fptr, "%d", &m);
    
    supMatSol = (double*)malloc(m*sizeof(double));
    if(supMatSol==NULL){
        printf("Memoria insuficient.\n");
        exit(1);
    }
    for(i=0; i<m; i++){
        supMatSol[i] = 0;
    }
    
    for(i=0; i<m; i++){
        for(j=0; j<n; j++){
            fscanf(fptr, "%lf", &termInd[j]);
        }
        fprintf(sort, "Vector b_%d (termes independents):\n", i+1);
        fprintf(sort, "(");
        for(j=0; j<n-1; j++){
            fprintf(sort, " %12.9e,", termInd[j]);
        }
        fprintf(sort, " %12.9e", termInd[n-1]);
        fprintf(sort, " )\n");
        
        /*
        permutar terminos independientes metodo SUSANA
        for(j=0; j<n; j++){
            if(perm[j] > j){
                aux = termInd[j];
                termInd[j] = termInd[perm[j]];
                termInd[perm[j]] = aux;
            }
        }
         */
        
        
        /*METODO MIO*/
        auxV = (double*)malloc(n*sizeof(double));
        if(auxV==NULL){
            printf("Memoria insuficient.\n");
            exit(1);
        }
        for(i=0; i<n; i++){
            auxV[i] = termInd[i];
        }
        for(i=0; i<n; i++){
            termInd[i] = auxV[perm[i]];
        }
        /*------*/
        
        triang(n, a, termInd, -1);
        triang(n, a, termInd, 1);
        
        fprintf(sort, "Vector x_%d (solucio):\n", i+1);
        fprintf(sort, "(");
        for(j=0; j<n-1; j++){
            fprintf(sort, " %12.9e,", termInd[j]);
        }
        fprintf(sort, " %12.9e", termInd[n-1]);
        fprintf(sort, " )\n\n");
        for(j=0; j<n; j++){
            supMatSol[j] += fabs(termInd[j]);
        }
    }
    
    for(i=0; i<m; i++){
        if(supMatSol[supInd] < supMatSol[i]){
            supInd = i;
        }
    }
    
    fprintf(sort, "S'han realitzat %d permutacions en la factoritzacio LU.\n", nP);
    fprintf(sort, "Norma del suprem de la matriu d'entrada = %12.9lf\n", supMatSist);
    fprintf(sort, "Norma del suprem de la matriu solucio = %12.9lf\n", supMatSol[supInd]);
    
    fclose(fptr);
    fclose(sort);
    
    for(i=0; i<n; i++){
        free(a[i]);
    }
    free(a);
    free(perm);
    
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

int lu(double **a, int n, int *perm, double tol){
    int i, j, k, maxIndex, nP=0, aux;
    double *auxV;
    
    for(j=0; j < (n-1); j++){
        /*buscar maximo*/
        maxIndex = j;
        /*buscar maximo de la columna*/
        for(i=j+1; i<n; i++){
            if( fabs(a[i][j])>fabs(a[maxIndex][j]) ){
                maxIndex = i;
            }
        }
        if( fabs(a[maxIndex][j]) < tol ){
            printf("Error: divisió per zero.");
            return -1;
        }
        
        if(maxIndex!=j){
            /*vector permutacion*/
            aux = perm[j];
            perm[j] = perm[maxIndex];
            perm[maxIndex] = aux;

            /*pivotar filas*/
            auxV = a[j];
            a[j] = a[maxIndex];
            a[maxIndex] = auxV;
            
            nP++;
        }
            
        /*multiplicadores*/
        for(i=j+1; i<n; i++){
            a[i][j] = a[i][j]/a[j][j];
        }
        
        /*hacer ceros*/
        for(i=j+1; i<n; i++){
            for(k=j+1; k<n; k++){
                a[i][k] = a[i][k]-a[i][j]*a[j][k];
            }
        }
        
    }
    
    return nP;
}
