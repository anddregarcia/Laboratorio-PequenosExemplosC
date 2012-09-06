//algoritmo Torre de hanoi

//arquivos de cabecalho
#include <iostream.h>
#include <conio.h>
#include <graphics.h>
#include <dos.h>
#include <stdio.h>
#include <stdlib.h>


int pilha[3][10];
int altura[3] = {0,0,0};
int disco[10] = {1,2,3,4,5,6,7,8,9,10};

// coloca o disco na torre
void empilha(int x, int torre)
{
   pilha[torre][altura[torre]] = x;
   altura[torre]++;
}

// remove o disco da torre
int desempilha(int torre)
{
    altura[torre]--;
    return pilha[torre][altura[torre]];
}

/*  parte grafica */
//inicia o grafico
int inicia_grafico()
{
    int gdriver = DETECT, gmode;
    initgraph(&gdriver, &gmode, "c:\\borlandc\\bgi");

    return graphresult();
}

//retira o disco graficamente
void apaga_disco(int tam,int torre)
{
    int x, y;

    y  = (getmaxy() / 2) - ((altura[torre] + 1) * 10);
    x  = (getmaxx() / 2) + ((torre - 1) * 100);

    setcolor(BLACK);

    tam = 10 - tam;
    tam *= 5;

    line(x - tam, y,x - tam,y + 9);
    line(x - tam, y,x + tam,y);
    line(x + tam, y,x + tam,y + 9);
}

//coloca o disco graficamente
void desenha_disco(int tam,int torre)
{
    int x, y;

    x  = (getmaxx() / 2) + ((torre - 1) * 100);
    y  = (getmaxy() / 2) - ((altura[torre]) * 10);

    tam = 10 - tam;
    tam *= 5;

    setcolor(WHITE);
    rectangle(x - tam,y,x + tam,y+10);
}

//algoritmo de movimentacao de disco
void movediscos(int n, int torreI, int torreF)
{
    //variaveis locais
    int torreAux = 3 - (torreI + torreF);
    int tam;

    //se existe apenas um disco a ser movido
    if(n == 1)
    {
     tam = desempilha(torreI);
     empilha(tam,torreF);

     apaga_disco(tam,torreI);
     desenha_disco(tam,torreF);

     delay(500);
    }
    else
    {
    movediscos(n-1,torreI, torreAux);
    movediscos(1, torreI, torreF);
    movediscos(n-1,torreAux,torreF);
    }
}

void main()
{
    int i, y, tam_torre;

    cout << "\nEntre com o tamanho da torre: ";
    cin  >> tam_torre;

    if(tam_torre > 10)
    {
     cout << "\n O Programa aceita no maximo 10 discos";
     exit(1);
    }

    inicia_grafico();
    line(1,getmaxy() / 2 ,getmaxx(),getmaxy() / 2);

    for(i=0;i<tam_torre;i++)
    {
     altura[0]++;
     pilha[0][i] = disco[i];
     desenha_disco(i + 1,0);
    }

    movediscos(tam_torre,0,1);

    getch();
}