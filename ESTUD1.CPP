#include <stdio.h>
#include <conio.h>

void main(void)
{
   struct
	 { int cod;
	   char nome[30];
			 } tabela[5];
   int i, coddesejado;


   clrscr();

   for(i=1;i<=5;i++)
   {
      printf("\n\n\tO c¢digo: ");
      scanf("%d",tabela[i].cod);
      flushall();
      printf("\n\n\tO nome: ");
      gets(tabela[i].nome);

     }





   getch();
}