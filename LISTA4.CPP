#include <stdio.h>
#include <conio.h>
#include <malloc.h>

struct reg {
  int id;
  struct reg *prox;
  struct reg *ant;
};

struct reg *pos = NULL;

void add(int novo_id) {
  struct reg *novo = (struct reg *) malloc(sizeof(reg));
  novo->id=novo_id;
  if (pos==NULL) {
     novo->prox=NULL; // A
     novo->ant=NULL;  // B
     pos=novo;        // C
  } else if (pos->prox!=NULL) {
    novo->prox=pos->prox;  // A
    novo->ant=pos;         // B
    pos->prox=novo;        // C
    novo->prox->ant=novo;  // D
  } else {
    novo->prox=pos;  // A
    novo->ant=pos;   // B
    pos->prox=novo;  // C
    pos->ant=novo;   // D
  }
}

void rmv() {
  if (pos!=NULL) {
    if (pos->prox!=NULL) {
      pos->ant->prox=pos->prox;
      pos->prox->ant=pos->ant;
      if (pos==pos->ant || pos==pos->prox) {
	free(pos);
	pos=NULL;
      } else {
	struct reg *rmv = pos;
	pos=pos->ant;
	free(rmv);
      }
    } else {
	free(pos);
	pos=NULL;
    }
  }
}

void lista() {
  struct reg *atual = pos;
  if (atual!=NULL) {
    do {
      printf("%d\n",atual->id);
      atual=atual->prox;
    } while (atual!=pos);
  } else {
    printf("LISTA VAZIA!");
  }
}

void main(void)
{
   clrscr();
   add(1);
   add(2);
   add(3);
   add(4);
   lista();
   printf("Pos = %d\n\n",pos->id);
   rmv();
   lista();
   printf("Pos = %d\n\n",pos->id);
   rmv();
   lista();
   printf("Pos = %d\n\n",pos->id);
   rmv();
   lista();
   printf("Pos = %d\n\n",pos->id);
   rmv();
   lista();
   getch();
}