

/*

Lepore factorization nr. 105 (Bruteforce) is free software
Lepore factorization nr. 105 (Bruteforce) is Alberico Lepore creation
this file is free software

The Lepore factorization nr. 105 (Bruteforce) is free software; you can redistribute it and/or modify
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
void MCD_Euclide(mpz_t A,mpz_t B , mpz_t *MCD);

int main(){

/*
This algorithm is generic and does not exploit that q / p < 2

Plus it uses a single A and not many A

A=9+16*a;//choose A with many small factors

default A =105 

to see "set this value"

*/




mpf_set_default_prec(1000000);/*set this value*/

mpz_t N,zero,uno,due,tre,quattro,cinque,otto,mom1,mom2,mom3,copia_N,p,q,A,M,x,P,Q,a;

mpf_t mom1_f,mom2_f,zero_f,uno_f,due_f,quattro_f,x_f,M_f,P_f,Q_f;
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
mpz_init(mom3);
mpz_init_set_str (zero, "0", 10);
mpz_init_set_str (uno, "1", 10);
mpz_init_set_str (due, "2", 10);
mpz_init_set_str (tre, "3", 10);
mpz_init_set_str (quattro, "4", 10);
mpz_init_set_str (cinque, "5", 10);
mpz_init_set_str (otto, "8", 10);
mpz_init(M);
mpz_init(x);
mpz_init(P);
mpz_init(Q);
mpz_init_set_str (a, "3", 10);
mpz_init_set_str (A, "105", 10);/*set this value*/

mpf_init(mom1_f);
mpf_init(mom2_f);
mpf_init_set_str (zero_f, "0.0", 10);
mpf_init_set_str (uno_f, "1.0", 10);
mpf_init_set_str (due_f, "2.0", 10);
mpf_init_set_str (quattro_f, "4.0", 10);
mpf_init(x_f);
mpf_init(M_f);
mpf_init(P_f);
mpf_init(Q_f);



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



while(1){
mpf_set_z(M_f,M);
mpf_set_z(mom1_f,M);
mpf_add(mom1_f,mom1_f,uno_f);
mpf_sqrt(mom1_f,mom1_f);
mpf_add(mom1_f,mom1_f,due_f);
mpf_div(mom1_f,mom1_f,quattro_f);

mpz_set_f(mom1,mom1_f);
mpf_set_z(mom2_f,mom1);

if(mpf_cmp(mom1_f,mom2_f)!=0){
mpf_set(x_f,mom2_f);
}else{
mpf_sub(x_f,mom2_f,uno_f);
}

mpf_mul(mom1_f,due_f,x_f);
mpf_add(mom1_f,mom1_f,uno_f);
mpf_mul(mom1_f,mom1_f,mom1_f);
mpf_mul(mom1_f,quattro_f,mom1_f);
mpf_sub(mom1_f,mom1_f,M_f);
if(mpf_cmp(mom1_f,zero_f)>=0){//if1
mpf_sqrt(mom1_f,mom1_f);
mpf_mul(mom2_f,quattro_f,x_f);
mpf_add(mom2_f,mom2_f,due_f);
mpf_sub(P_f,mom2_f,mom1_f);

mpz_set_f(P,P_f);
mpf_set_z(mom2_f,P);



if(mpf_cmp(P_f,mom2_f)==0){//if2
mpz_div(Q,M,P);

mpz_mod(mom1,P,N);
mpz_mod(mom2,Q,N);

if(mpz_cmp(mom1,zero)!=0 && mpz_cmp(mom2,zero)!=0){
break;
}

}//if2
}//if1


mpz_mul(M,M,A);

mpz_add(a,a,quattro);
}//while(1)


MCD_Euclide(N,P,&p);
gmp_printf ("\np=%Zd\n",p);

}



void MCD_Euclide(mpz_t A,mpz_t B , mpz_t *MCD){

mpz_t R,zero,mom1;

mpz_init(R);
mpz_init_set_str (zero, "0", 10);
mpz_init(mom1);

while(mpz_cmp(B,zero)!=0){

if(mpz_cmp(B,A)>0){
mpz_set(mom1,A);
mpz_set(A,B);
mpz_set(B,mom1);
}

mpz_mod(R,A,B);
mpz_set(A,B);
mpz_set(B,R);

}

mpz_set(*MCD,A);
return;
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
















