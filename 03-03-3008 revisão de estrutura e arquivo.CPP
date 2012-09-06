#include <stdio.h>
#include <conio.h>

struct dado
{
  int id;
  char nome[30];
};

int abrir(struct dado info[5]) {
  FILE *arquivo;

  if ((arquivo = fopen("C:\\TESTE.DAT", "rb")) == NULL)
{
   fprintf(stderr, "Erro ao abrir o arquivo.\n");
   return 0;
}
fseek(arquivo, SEEK_SET, 0);
fread(&info, sizeof(dado), 5, arquivo);
fclose(arquivo);
return 1;

}

int salvar(struct dado info[5] ) {
   FILE *arquivo;
   if ((arquivo = fopen("C:\\TESTE.DAT", "wb")) == NULL)
   {
      fprintf(stderr, "Erro ao abrir o arquivo.\n");
      return 0;
   }
   fwrite(&info, sizeof(dado), 5, arquivo);
   fclose(arquivo);
   return 1;
}

void inserir(struct dado info[5]) {
   clrscr();
   for (int i=0; i<5; i++) {   
     printf("%d. Id: ",i); scanf("%d%c",&info[i].id);
     printf("%d. Nome: ",i); gets(info[i].nome);
     printf("\n");
   }
   getch();
}

void exibir(struct dado info[5]) {
   clrscr();
   for (int i=0; i<5; i++) {
     printf("%d. Id: %d ",i,info[i].id);
     printf(" Nome: %s\n",info[i].nome);
   }
   getch();
}


void main(void)
{
  struct dado info[5];
  char op;
  do {
    clrscr();
    printf("MENU PRINCIPAL\n");
    printf("--------------\n");
    printf("1. Inserir dados\n");
    printf("2. Visualizar dados\n");
    printf("3. Ler arquivo\n");
    printf("4. Gravar arquivo\n");
    printf("5. Sair\n");
    op=getch();
    if (op=='1') {
      inserir(info);
    } else if (op=='2') {
      exibir(info);
    } else if (op=='3') {
      clrscr();
      if (abrir(info)) {
	printf("Arquivo aberto com sucesso!");
      }
      getch();
    } else if (op=='4') {
      clrscr();
      if (salvar(info)) {
	printf("Arquivo salvo com sucesso!");
      }
      getch();
    }
  } while (op!='5');

}

