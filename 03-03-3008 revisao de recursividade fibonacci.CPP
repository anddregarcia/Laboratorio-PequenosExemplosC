#include <stdio.h>
#include <conio.h>

int fib(int f) {
  if (f<=1) {
    return f;
  } else {
    return fib(f-2) + fib(f-1);
  }
}

void main() {
  clrscr();
  printf("Fibonacci: %d",fib(9));
  getch();
}