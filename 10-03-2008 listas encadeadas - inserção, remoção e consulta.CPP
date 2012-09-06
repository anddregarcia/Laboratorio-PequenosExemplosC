#include <stdio.h>
#include <conio.h>
#include <stdlib.h>

struct reg {
	   int pid, pri, tipo, cpu, io;

  struct reg *prox;
};

struct esc {
	   
	   int pronto, exe, bloq;
	   
	   }
	   

struct reg *inicio = NULL;
struct reg *ultimo = NULL;

void adiciona(int pid, int pri, int tipo, int cpu, int io) {

  struct reg *novo;

  novo=(reg *) malloc(sizeof(reg));

  if (inicio==NULL) {

    inicio=novo;

  } else {

    ultimo->prox=novo;

  }

  novo->pid=pid;
  novo->pri=pri;
  novo->tipo=tipo;
  novo->cpu=cpu;
  novo->io=io;

  novo->prox=NULL;

  ultimo=novo;

}

void apaga() {
  struct reg *atual;
  atual=inicio;
  if (inicio==ultimo) {
    free(ultimo);
    inicio=NULL;
    ultimo=NULL;
  } else {
    while (atual->prox!=ultimo) {
      atual=atual->prox;
    }
    free(ultimo);
    ultimo=atual;
    ultimo->prox=NULL;
  }
}

void limpar() {
  while (inicio!=NULL) {
    apaga();
  }
}

void lista() {
  struct reg *atual;
  atual=inicio;
  while (atual!=NULL) {
    printf("%d\t",atual->pid);
    printf("%d\t",atual->pri);
    printf("%d\t",atual->tipo);
    printf("%d\t",atual->cpu);
    printf("%d\n",atual->io);
    
    atual=atual->prox;
  }
}

void roundrobin(int pid, int cpu, int io){
	 struct reg *atual;
	 atual=inicio;
	 while (atual != NULL){
	 	   
	 	   if(atual->pid < pid->prox)
	 	   
	 
	 
	 







main() {
  int op;
  do {


    printf("MENU - ESCALONADOR DE PROCESSOS\n\n");
    printf("----------------------\n");
    printf("Processos-----\n");
    printf("1. Adicionar\n");
    printf("2. Remover\n");
    printf("3. Listar\n");
    printf("4. Limpar\n\n");
    printf("Escalonadores ------\n");
    printf("5. Round Robin\n\n\n");


    printf("0. Sair\n\n");
    printf("Opcao: "); scanf("%d",&op);

    switch (op) {
      case 1: {
	  
	  int pid, pri, tipo, cpu, io;

	printf("Digite o PID: ");
	scanf("%d",&pid);
	printf("Digite o PRI: ");
	scanf("%d",&pri);

	while((tipo>3)||(tipo<1))
	{
	printf("Digite o TIPO (1-system/2-user/3-batch): ");
	scanf("%d",&tipo);
	}

	printf("Digite o CPU: ");
	scanf("%d",&cpu);
	printf("Digite o I/O: ");
	scanf("%d",&io);

	adiciona(pid, pri, tipo, cpu, io);


	break;
      }
      case 2: {
	if (inicio!=NULL) {
	  printf("Removendo... ");
	  apaga();
	  getch();
	}
	break;
      }
      case 3: {
	lista();
        getch();
	break;
      }
      case 4: {
	printf("Limpando...");
	limpar();
	getch();
	break;
      }
    
    case 5: {
		 
		 roundrobin(pid, cpu, io);
		 
		 getch();
		 break;
		 }
		 
    
	}
    
  } while (op!=0);
  limpar();
}
