// Listas Encadeadas - Exemplo 1

#include <stdio.h>
#include <conio.h>
#include <alloc.h>

#define numreg 10

struct reg {
  int info;
  struct reg *prox;
};

struct reg *inicio = NULL;
struct reg *ultimo = NULL;

void adiciona(int info) {
  struct reg *novo;
  novo=(reg *) malloc(sizeof(reg));
  if (inicio==NULL) {
    inicio=novo;
  } else {
    ultimo->prox=novo;
  }
  novo->info=info;
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
    printf("%d\n",atual->info);
    atual=atual->prox;
  }
}

void main() {
  int op;
  do {
    clrscr();
    printf("MENU - LISTA ENCADEADA\n");
    printf("----------------------\n");
    printf("1. Adicionar\n");
    printf("2. Remover\n");
    printf("3. Listar\n");
    printf("4. Limpar\n");
    printf("5. Sair\n\n");
    printf("Opcao: "); scanf("%d",&op);
    clrscr();
    switch (op) {
      case 1: {
	printf("Digite o dado: ");
	int dado;
	scanf("%d",&dado);
	adiciona(dado);
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
    }
  } while (op!=5);
  limpar();
}