#include <cstdio>
#include <cstdint>
#include <vector>
#include "ligra.h"

int main(int argc, char **argv){
    uint64_t x;
    commandLine P(argc, argv, " [-s] <inFile>");
    if (P.getOption("-n"))
    {
        uint64_t num_vertices = P.getOptionLongValue("-n");
        fwrite(&num_vertices, sizeof(num_vertices), 1, stdout);
    }
    while(scanf("%lu", &x) != EOF){
        fwrite(&x, sizeof(x), 1, stdout);
    }
}
