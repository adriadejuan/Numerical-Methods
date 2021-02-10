#include<stdio.h>
#include<math.h>

double eTaylor0(int, double);

int main(void){
    int i = 0, grau, fact = 1;
    double x, sol, errorFitat, errorMaxim = (1e-8)/2;
    
    printf("X = 1\nX_0 = 0\nGrau: ");
    scanf("%d",&grau);
    printf("\n");
    printf("               Grau:              Valor:        Error fitat:\n");
    
    x = 1;
    
    sol = eTaylor0(i, x);
    errorFitat = 2.8;
    printf("%20.d%20.10e%20.10e\n",i,sol,errorFitat);
    
    for (i = 1; i<=grau; i++) {
        sol = eTaylor0(i, x);
        errorFitat = fabs(2.8/fact);
        printf("%20.d%20.10e%20.10e\n",i,sol,errorFitat);
        fact = fact * i;
    }
    
    i = 1;
    fact = 1;
    while ((2.8/fact)>errorMaxim) {
        i++;
        fact = fact * i;
    }
    
    printf("\nEs poden garantir 8 xifres decimals correctes amb n = %d\n",i);
    
    return 0;
}

double eTaylor0(int grau, double x){
    int i, fact=1;
    double prod, ans;
    
    ans = 1;
    
    if (grau==0){
        return ans;
    }
    
    prod = x;
    
    for(i=1;i<=grau;i++){
        fact = fact * i;
        ans += prod/fact;
        prod = prod*x;
    }
    
    return ans;
}
