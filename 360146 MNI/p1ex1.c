#include<stdio.h>
#include<math.h>

int main(void){
    int i, k, iter;
    double s, sAnt, c, error = 999;
    
    k = 4;
    s = sqrt(2);
    c = 2*s/sqrt(4-s*s);
    
    printf("Iteracions: ");
    scanf("%d",&iter);
    
    printf("         Costats(n):       Inscrit(S_k):   Circumscrit(C_k):             e_r(π):\n");
	
	for(i=1; i<=iter; i++){
        printf("%20.d%20.8e%20.8e%20.8e\n",k,s*k/2,c*k/2,error);
        k = k*2;
        sAnt = s;
        s = sqrt(2-sqrt(4-s*s));
        error = (s - sAnt)/s;
        c = 2*s/sqrt(4-s*s);
	}
    
    return 0;
}
