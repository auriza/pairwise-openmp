//  columnwise per n/p block; manual; loop unroll 2x

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

//#define DEBUG
#define MATCH       +1
#define MISMATCH    -1
#define GAP         -3

double time;

void LCS_length(const char *X, const char *Y, short int ***c_out, char ***b_out)
{
    int m;          // row
    int n;          // column
    short int **c;  // length matrix
    char **b;       // arrow matrix
    int diag, left, up; // diagonal (upper-left), left, up score
    int i, j;

    m = strlen(X);
    n = strlen(Y);
    c = malloc((m+1)*sizeof(short int *));
    b = malloc((m+1)*sizeof(char *));
    for (i = 0; i <= m; i++) {
        c[i] = malloc((n+1)*sizeof(short int));
        b[i] = calloc(n+1, sizeof(char));
    }

    // init
    c[0][0] = 0;
    b[0][0] = ' ';
    for (i = 1; i <= m; i++) {
        c[i][0] = i * GAP;
        b[i][0] = '|';
    }
    for (j = 1; j <= n; j++) {
        c[0][j] = j * GAP;
        b[0][j] = '_';
    }

    time = omp_get_wtime();

    // fill the matrix
    #pragma omp parallel private(i,j,diag,left,up)
    {
        int id = omp_get_thread_num();
        int t  = omp_get_num_threads();
        int jfirst = id*n/t + 1;
        int jlast  = (id+1)*n/t;

        for (i = 1; i <= m; i += 2) {

            // check if leftmost cell value is not empty?
            while (b[i][jfirst-1] == 0);

            for (j = jfirst; j <= jlast; j++) {
                // calculate scoring possibility
                diag = c[i-1][j-1] + ((X[i-1] == Y[j-1])? MATCH : MISMATCH);
                up   = c[i-1][j] + GAP;
                left = c[i][j-1] + GAP;

                if (diag >= up && diag >= left) {
                    c[i][j] = diag;
                    b[i][j] = '\\';
                } else if (up >= left) {
                    c[i][j] = up;
                    b[i][j] = '|';
                } else {
                    c[i][j] = left;
                    b[i][j] = '_';
                }

                if (i == m) break;

                diag = c[i][j-1] + ((X[i] == Y[j-1])? MATCH : MISMATCH);
                up   = c[i][j] + GAP;
                left = c[i+1][j-1] + GAP;

                if (diag >= up && diag >= left) {
                    c[i+1][j] = diag;
                    b[i+1][j] = '\\';
                } else if (up >= left) {
                    c[i+1][j] = up;
                    b[i+1][j] = '|';
                } else {
                    c[i+1][j] = left;
                    b[i+1][j] = '_';
                }
            }
        }
    }

    time = omp_get_wtime() - time;

    *c_out = c;
    *b_out = b;
}

void LCS_print_X(char **b, const char *X, int i, int j)
{
    if (i == 0 && j == 0) {
        return;
    }
    if (b[i][j] == '\\') {
        LCS_print_X(b, X, i-1, j-1);
        printf("%c", X[i-1]);
    } else if (b[i][j] == '|') {
        LCS_print_X(b, X, i-1, j);
        printf("%c", X[i-1]);
    } else {
        LCS_print_X(b, X, i, j-1);
        printf("_");

    }
}

void LCS_print_Y(char **b, const char *Y, int i, int j)
{
    if (i == 0 && j == 0) {
        return;
    }
    if (b[i][j] == '\\') {
        LCS_print_Y(b, Y, i-1, j-1);
        printf("%c", Y[j-1]);
    } else if (b[i][j] == '|') {
        LCS_print_Y(b, Y, i-1, j);
        printf("_");
    } else {
        LCS_print_Y(b, Y, i, j-1);
        printf("%c", Y[j-1]);
    }
}

void LCS_print_c(short int **c, int m, int n)
{
    int i, j;
    for (i = 0; i <= m; i++) {
        for (j = 0; j <= n; j++)
            printf("%3hd ", c[i][j]);
        printf("\n");
    }
    printf("\n");
}

void LCS_print_b(const char X[], const char Y[], char **b, int m, int n)
{
    int i, j;
    printf(" ");
    for (i = 1; i <= n; i++)
        printf(" %c", Y[i-1]);
    printf("\n");
    for (i = 1; i <= m; i++) {
        printf("%c", X[i-1]);
        for (j = 1; j <= n; j++) {
            printf(" %c", b[i][j]);
        }
        printf("\n");
    }
}

int main(void)
{
    short int **c;  // score matrix
    char **b;       // direction matrix
    char *X, *Y;    // sequences
    int m, n;       // sequences' length

    scanf("%d", &m);
    scanf("%d", &n);
    X = malloc((m+1) * sizeof(char));
    Y = malloc((n+1) * sizeof(char));
    scanf("%s", X);
    scanf("%s", Y);

    LCS_length(X, Y, &c, &b);

    #ifdef DEBUG
        // print LCS
        LCS_print_X(b, X, m, n); printf("\n");
        LCS_print_Y(b, Y, m, n); printf("\n\n");

        LCS_print_c(c, m, n);
        LCS_print_b(X, Y, b, m, n);
    #endif

    printf("Time: %f s\tScore: %d\n", time, c[m][n]);

    return 0;
}