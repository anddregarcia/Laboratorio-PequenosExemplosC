#include <stdlib.h>
#include <stdio.h>
#include <conio.h>

int main(void)
{
  float preco[20];
  float mpn,mpd;
  int i,fpn,fpd;

  i=0;
  while (i<19)
  {
    printf ("Fornecedor %d",i/2+1);
    printf ("Preco Normal: ");
    scanf ("%f",&preco[i]);
    printf ("Preco Desconto: ");
    scanf ("%f",&preco[i+1]);
    i+=2;
  }
  mpn=preco[0];
  fpn=1;
  mpd=preco[1];
  fpd=1;
  i=0;
  while (i<19)
  {
    if (preco[i]>mpn)
    {
      mpn=preco[i];
      fpn=i;
    }
    if (preco[i+1]>mpd)
    {
      mpd=preco[i+1];
      fpd=1;
    }
    i=i+2;
  }

  printf ("Maior Preco Normal %5.2f",mpn);
  printf("  Fornecedor: %d",fpn/2+1);

  printf ("Maior Preco Desconto %5.2f",mpd);
  printf ("  Fornecedor: %d",fpd/2+1);
  getch();
}
