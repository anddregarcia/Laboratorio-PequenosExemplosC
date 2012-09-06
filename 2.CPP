// 0504 - ANDRê DE SOUZA ANANIAS GARCIA

#include <stdio.h>
#include <conio.h>

void main()
{  struct { int salario; int idade; int filhos; char sexo; } vet[100];

   int i;
   int sal_tot, media_sal;
   int filho_tot, media_filho;
   int porcent, sal_mu;
   int cont;
   int habit;

 clrscr();


   printf("\n\n\tDigite a quantidade de habitantes pesquisados: ");
   scanf("%d",&habit);

   for(i=0;i<=habit;i++)
   {
    clrscr();

    printf("\n\t\t---LEITURA DOS DADOS DOS HABITANTE---");

    printf("\n\n\n\tDigite o sal†rio: R$ ");
    scanf("%d",&vet[i].salario);

    while (vet[i].salario<0)
    {
      printf("\n\n\tSALARIO INVALIDO!!!");
      printf("\n\n\tDigite novamente: R$ ");
      scanf("%d",&vet[i].salario);
      clrscr();
    }

    flushall();
    printf("\n\n\tDigite o sexo (M/F): ");
    scanf("%c",&vet[i].sexo);

    printf("\n\n\tDigite a idade: ");
    scanf("%d",&vet[i].idade);

    printf("\n\n\tDigite a quantidade de filhos: ");
    scanf("%d",&vet[i].filhos);

   }
    clrscr();

   printf("\t\tDADOS CALCULADOS");

//MEDIA DO SALARIO

   sal_tot=0;

   for(i=0;i<=habit;i++)
   {
      sal_tot=(sal_tot + vet[i].salario);

      }

      media_sal=(sal_tot/habit);
      printf("\n\n\n\tA mÇdia do salario dos habitantes Ç: R$ %.2d", media_sal);

//MEDIA DE Nß DE FILHOS

   filho_tot=0;

   for(i=0;i<=habit;i++)
   {
      filho_tot=(filho_tot + vet[i].filhos);

     }

      media_filho=(filho_tot/habit);
      printf("\n\n\n\tA mÇdia da quantidade de filhos de todo populaáao Ç: %.2d", media_filho);

//PERCENTUAL DAS MELHERES COM SALARIO MAIOR QUE R$ 1.000,00

    cont=0;

    for(i=0;i<=habit;i++)
    {

      if((vet[i].sexo=='f')||(vet[i].sexo=='F'))
      {
	if(vet[i].salario>1000)
	{
	   cont=cont+1;

	    }
	 }

      }

       printf("\n\n\n\tH† %d mulheres na cidade.",cont);
       porcent=(100*cont/habit);
       printf("\n\n\n\tO percentual das mulheres que ganham mais que R$1.000,00 Ç: %.2d % ",porcent);


   getch();

}
