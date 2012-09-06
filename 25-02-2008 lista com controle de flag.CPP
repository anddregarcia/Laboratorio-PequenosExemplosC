#include <stdio.h>
#include <conio.h>

int lista[10][2];

void set(int pos) {
  lista[pos][0]=0;
}

void out(int pos) {
  lista[pos][0]=-1;
}

void preenche() {
  clrscr();
  for (int i=0; i<10; i++) {
    printf("[%d]: ",i);
    int tmp;
    scanf("%d",&lista[i][1]);
  }
}

void list() {
  clrscr();
  for (int i=0; i<10; i++) {
    if (lista[i][0]==0) {
      printf("[%d] %d\n",i,lista[i][1]);
    }
  }
  getch();
}


void main(void)
{
   int opcao;
   do {
     clrscr();
     printf("MENU\n");
     printf("------\n");
     printf("(1) Preenche\n");
     printf("(2) Set-In\n");
     printf("(3) Set-Out\n");
     printf("(4) Lista\n");
     printf("(5) Sair\n");
     printf("\nOpcao: ");
     scanf("%d",&opcao);
     switch (opcao) {
       case 1: {
	 preenche();
	 break;
       }
       case 2: {
	 int id;
	 do {
	   clrscr();
	   printf("Set-In\nIdentificacao (0 a 9): ");
	   scanf("%d",&id);
	 } while (!(id>=0 && id<=9));
	 set(id);
	 break;
       }
       case 3: {
	 int id;
	 do {
	   clrscr();
	   printf("Set-Out\nIdentificacao (0 a 9): ");
	   scanf("%d",&id);
	 } while (!(id>=0 && id<=9));
	 out(id);
	 break;
       }
       case 4: {
	 list();
         break;
       }
     }
   } while (opcao!=5);
   getch();
}