#include <stdio.h>
#include <conio.h>

int fatorial(int f) {
  if (f==0) { return 1; }
  else { return f * fatorial(f-1); }
}

void main() {
  clrscr();
  printf("Fatorial: %d",fatorial(5));
  getch();
}