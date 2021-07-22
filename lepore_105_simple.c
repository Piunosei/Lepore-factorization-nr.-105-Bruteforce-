/*
Lepore factorization nr. 105 simple  is free software
Lepore factorization nr. 105 simple is Alberico Lepore creation
this file is free software
The Lepore factorization nr. 105 simple is free software; you can redistribute it and/or modify
it under the terms of either:
    the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
You should have received copies of the
GNU Lesser General Public License here https://www.gnu.org/licenses/.
*/







/* gmp_version -- 6.1.2
Copyright 1996, 1999-2001 Free Software Foundation, Inc.
This file is part of the GNU MP Library.
The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of either:
  * the GNU Lesser General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
or
  * the GNU General Public License as published by the Free Software
    Foundation; either version 2 of the License, or (at your option) any
    later version.
or both in parallel, as here.
The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.
You should have received copies of the GNU General Public License and the
GNU Lesser General Public License along with the GNU MP Library.  If not,
see https://www.gnu.org/licenses/.  */


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>


int prendi_numero(char in[]);


int main(){
/*
this algorithm factorizes 
N in the form  N=8*G+3=p*q
with
1 <= (q-p+2)/4 < (sqrt(32*(p+q-4)/8+1)+1)/2

*/
mpf_set_default_prec(1000);/*set this value*/

mpz_t N,copia_N,M,pseudo_M,m,zero,uno,due,tre,quattro,cinque,otto,p,q,x,y,mom1,mom2;

mpf_t mom1_f,m_f,x_f,due_f,uno_f;
mpz_init(N);

char numero_N[10000];/*set this value*/
prendi_numero(numero_N);   
mpz_init_set_str(N, numero_N, 10);

gmp_printf ("\nN=\n%Zd\n",N);


mpz_init(q);
mpz_init(p);
mpz_init(copia_N);
mpz_init(mom1);
mpz_init(mom2);
mpz_init_set_str (zero, "0", 10);
mpz_init_set_str (uno, "1", 10);
mpz_init_set_str (due, "2", 10);
mpz_init_set_str (tre, "3", 10);
mpz_init_set_str (quattro, "4", 10);
mpz_init_set_str (cinque, "5", 10);
mpz_init_set_str (otto, "8", 10);
mpz_init(M);
mpz_init(x);
mpz_init(y);
mpz_init(pseudo_M);
mpz_init(m);
mpz_init(mom1);
mpz_init(mom2);




mpf_init(mom1_f);
mpf_init_set_str (uno_f, "1.0", 10);
mpf_init_set_str (due_f, "2.0", 10);
mpf_init(x_f);
mpf_init(m_f);


int check=0;


mpz_set(copia_N,N);
mpz_mod(mom1,N,quattro);
if(mpz_cmp(mom1,uno)==0){
mpz_mul(copia_N,copia_N,tre);
}

mpz_mod(mom1,copia_N,otto);
if(mpz_cmp(mom1,tre)!=0){
mpz_mul(M,copia_N,cinque);
}else{
mpz_set(M,copia_N);
}

mpz_sub(mom1,M,tre);
mpz_div(m,mom1,otto);

mpf_set_z(m_f,m);
mpf_mul(mom1_f,m_f,due_f);
mpf_add(mom1_f,mom1_f,uno_f);
mpf_sqrt(mom1_f,mom1_f);

mpf_sub(mom1_f,mom1_f,uno_f);
mpf_div(mom1_f,mom1_f,due_f);

mpz_set_f(x,mom1_f);
mpf_set_z(x_f,x);

if(mpf_cmp(x_f,mom1_f)!=0){
mpz_add(x,x,uno);
}

mpz_mul(mom1,x,x);
mpz_add(mom1,mom1,x);
mpz_mul(mom1,mom1,due);
mpz_sub(mom1,mom1,m);
mpz_mul(mom1,mom1,otto);
mpz_add(mom1,mom1,uno);
mpz_sqrt(mom1,mom1);
mpz_add(mom1,mom1,uno);
mpz_div(y,mom1,due);

mpz_mul(mom1,x,quattro);
mpz_add(mom1,mom1,tre);
mpz_mul(mom2,y,due);
mpz_sub(p,mom1,mom2);
mpz_div(q,M,p);
mpz_mul(pseudo_M,p,q);

if(mpz_cmp(p,uno)!=0 && mpz_cmp(p,M)!=0 && mpz_cmp(pseudo_M,M)==0){
check=1;;
}

if(check==1){
mpz_mod(mom1,p,tre);
if(mpz_cmp(mom1,zero)==0){
mpz_div(p,p,tre);
}

mpz_mod(mom1,p,cinque);
if(mpz_cmp(mom1,zero)==0){
mpz_div(p,p,cinque);
}
mpz_div(q,N,p);
gmp_printf ("\np=%Zd\nq=%Zd\n",p,q);
}else{
printf("\nN does not have these characteristics: N in the form  N=8*G+3=p*q with 1 <= (q-p+2)/4 < (sqrt(32*(p+q-4)/8+1)+1)/2\n");
}


}






int prendi_numero(char in[]){

    char decimale[1000];
    int numero_di_cifre_decimali=0;
    FILE *fp;
    int i=0;

    fp = fopen("input.txt", "r");
    if (fp==NULL){
        printf("\nImpossibile aprire file\n");
        system("PAUSE");
        exit(1);
    }
    while(!feof(fp)){
        fscanf(fp,"%s",decimale);

    }
    fclose(fp);

    numero_di_cifre_decimali=strlen(decimale)-1;
    while(i<=numero_di_cifre_decimali){
        in[i]=decimale[i];
        i++;
    }
    return numero_di_cifre_decimali;
}
