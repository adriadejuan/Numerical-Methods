#include<stdio.h>
#include<math.h>

int main(void){
	int i, n, signe=-1;
	double iAnt, frac, e0, eN;
    iAnt = 0.125657;
    frac = (7./2);
    e0 = -0.00000021;
    eN = signe * frac * e0;
	
	printf("Entra n: ");
	scanf("%d",&n);
	
    printf("                  n:                I_n:                e_n:\n");
    printf("%20.d%20.10e%20.10e\n",0,iAnt,e0);
    
	for(i=1; i<=n; i++){
        eN = signe * frac * e0;
        printf("%20.d%20.10e%20.10e\n",i,iAnt,eN);
        signe = signe * (-1);
        frac = frac * (7./2);
	}
    return 0;
}
