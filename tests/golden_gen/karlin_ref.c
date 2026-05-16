/*
 * Standalone reference harness for Karlin/Altschul window-size estimation.
 *
 * Extracts karlin() and winsizeFromlambdak() from
 * third_party/seqtools/dotterApp/dotterKarlin.c verbatim (modulo a tiny
 * glib shim — we use plain malloc/fprintf instead of g_malloc/g_critical).
 * The numerical body is unchanged.
 *
 * Compile:
 *   gcc -O2 -o karlin_ref tests/golden_gen/karlin_ref.c -lm
 *
 * Output: tab-separated rows, one per built-in fixture, with columns
 *   name  lambda  k  h  exp_res_score  exp_msp_score  predicted  window
 *
 * The Rust port's goldens are pinned against this output.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int32_t gint32;

#define MAXIT 20
#define SUMLIMIT 0.01
#define NR 23
#define NA 24

/* ---- Helpers extracted from dotterKarlin.c verbatim ---- */

static double fct_powi(double x, int n) {
    int i;
    double y;
    y = 1.;
    for (i = abs(n); i > 0; i /= 2) {
        if (i & 1) y *= x;
        x *= x;
    }
    return n >= 0 ? y : 1. / y;
}

static long fct_gcd(long a, long b) {
    long c;
    b = labs(b);
    if (b > a) { c = a; a = b; b = c; }
    while (b != 0) { c = a % b; a = b; b = c; }
    return a;
}

static double fct_expm1(double x) {
    double absx = (x < 0) ? -x : x;
    if (absx > .33) return exp(x) - 1.;
    if (absx < 1.e-16) return x;
    return x * (1. + x *
        (0.5 + x * (1./6. + x *
            (1./24. + x * (1./120. + x *
                (1./720. + x * (1./5040. + x *
                    (1./40320. + x * (1./362880. + x *
                        (1./3628800. + x * (1./39916800. + x *
                            (1./479001600. + x/6227020800.)
                            ))
                        ))
                    ))
                ))
            )));
}

static double etop(double E) { return -fct_expm1(-E); }

/* ---- karlin(): direct port from dotterKarlin.c:144 ---- */

static double karlin(long low, long high, double *pr,
                     double *lambda, double *K, double *H) {
    int i, j;
    long range, lo, hi, first, last;
    double up, new_val, sum, Sum, av, beta, oldsum, oldsum2;
    double *p = NULL, *P = NULL, *ptrP, *ptr1, *ptr2;
    double ratio;

    if (low >= 0.) { fprintf(stderr, "karlin: low>=0\n"); return -1.0; }

    for (i = range = high - low; i > -low && pr[i] == 0.0; --i);
    if (i <= -low) { fprintf(stderr, "karlin: no positive score\n"); return -1.0; }

    for (sum = i = 0; i <= range; sum += pr[i++])
        if (pr[i] < 0.) { fprintf(stderr, "karlin: negative prob\n"); return -1.0; }

    p = (double *)malloc(sizeof(*p) * (range + 1));
    for (Sum = low, i = 0; i <= range; ++i)
        Sum += i * (p[i] = pr[i] / sum);

    if (Sum >= 0.) { fprintf(stderr, "karlin: non-neg expected score\n"); free(p); return Sum; }

    up = 0.5;
    do {
        up *= 2;
        ptr1 = p;
        for (sum = 0, i = low; i <= high; ++i)
            sum += *ptr1++ * exp(up * i);
    } while (sum < 1.0);

    for (*lambda = 0., j = 0; j < 25; ++j) {
        new_val = (*lambda + up) / 2.0;
        ptr1 = p;
        for (sum = 0., i = low; i <= high; ++i)
            sum += *ptr1++ * exp(new_val * i);
        if (sum > 1.0) up = new_val;
        else *lambda = new_val;
    }
    beta = exp(*lambda);

    ptr1 = p;
    for (av = 0, i = low; i <= high; ++i)
        av += *ptr1++ * i * exp(*lambda * i);
    *H = *lambda * av;

    if (low == -1 || high == 1) {
        *K = (high == 1 ? av : Sum * Sum / av);
        *K *= 1.0 - 1. / beta;
        free(p);
        return 0;
    }

    Sum = 0.;
    lo = hi = 0;
    P = (double *)malloc(MAXIT * (range + 1) * sizeof(*P));
    *P = sum = oldsum = oldsum2 = 1.;
    for (j = 0; j < MAXIT && sum > SUMLIMIT; oldsum = sum, Sum += sum /= ++j) {
        first = last = range;
        for (ptrP = P + (hi += high) - (lo += low); ptrP >= P; *ptrP-- = sum) {
            ptr1 = ptrP - first;
            ptr2 = p + first;
            for (sum = 0., i = first; i <= last; ++i)
                sum += *ptr1-- * *ptr2++;
            if (first != 0) --first;
            if (ptrP - P <= range) --last;
        }
        new_val = fct_powi(beta, lo - 1);
        for (sum = 0, i = lo; i != 0; ++i)
            sum += *++ptrP * (new_val *= beta);
        for (; i <= hi; ++i)
            sum += *++ptrP;
        oldsum2 = oldsum;
    }

    ratio = oldsum / oldsum2;
    if (ratio >= (1.0 - SUMLIMIT * 0.001)) {
        *K = 0.1;
        free(p); free(P);
        return 0;
    }
    while (sum > SUMLIMIT * 0.01) {
        oldsum *= ratio;
        Sum += sum = oldsum / ++j;
    }

    for (i = low; p[i - low] == 0.; ++i);
    for (j = -i; i < high && j > 1;)
        if (p[++i - low])
            j = fct_gcd(j, i);

    *K = (j * exp(-2. * Sum)) / (av * etop(*lambda * j));

    free(p); free(P);
    return 0;
}

/* ---- winsizeFromlambdak(): direct port from dotterKarlin.c:343 ---- */

static int winsizeFromlambdak(gint32 mtx[24][24], int *tob, int abetsize,
                              const char *qseq, const char *sseq,
                              double *exp_res_score, double *Lambda,
                              double *out_K, double *out_H, double *out_exp_msp) {
    gint32 lows = 0, highs = 0, range;
    int i, j;
    int *n1, *n2;
    int qlen = 0, slen = 0;
    int retval;
    int n = 100;
    double *fq1, *fq2, *prob, K, H;
    double qij, exp_MSP_score, sum;

    n1 = (int *)malloc((abetsize + 4) * sizeof(int));
    n2 = (int *)malloc((abetsize + 4) * sizeof(int));
    fq1 = (double *)malloc((abetsize + 4) * sizeof(double));
    fq2 = (double *)malloc((abetsize + 4) * sizeof(double));

    for (i = 0; i < abetsize; ++i)
        for (j = 0; j < abetsize; ++j) {
            if (mtx[i][j] < lows) lows = mtx[i][j];
            if (mtx[i][j] > highs) highs = mtx[i][j];
        }

    for (i = 0; i < abetsize; ++i) n1[i] = 0;
    for (i = 0; qseq[i]; ++i)
        if (tob[(int)(unsigned char)qseq[i]] < abetsize) {
            n1[tob[(int)(unsigned char)qseq[i]]]++;
            qlen++;
        }
    for (i = 0; i < abetsize; ++i) n2[i] = 0;
    for (i = 0; sseq[i]; ++i)
        if (tob[(int)(unsigned char)sseq[i]] != NA &&
            tob[(int)(unsigned char)sseq[i]] < abetsize) {
            /* The C original tests `!= NA` (=24); for DNA mode (abetsize=4)
             * this lets 'N' through to be filtered later. We also gate on
             * `< abetsize` to mirror the qseq loop. */
            n2[tob[(int)(unsigned char)sseq[i]]]++;
            slen++;
        }

    for (i = 0; i < abetsize; ++i) {
        fq1[i] = (double)n1[i] / qlen;
        fq2[i] = (double)n2[i] / slen;
    }

    range = highs - lows;
    prob = (double *)malloc(sizeof(double) * (range + 1));
    for (i = 0; i <= range; ++i) prob[i] = 0.0;
    for (i = 0; i < abetsize; ++i)
        for (j = 0; j < abetsize; ++j)
            prob[mtx[i][j] - lows] += fq1[i] * fq2[j];

    if ((*exp_res_score = karlin(lows, highs, prob, Lambda, &K, &H))) {
        free(prob); free(n1); free(n2); free(fq1); free(fq2);
        return 25;
    }

    *exp_res_score = sum = 0;
    for (i = 0; i < abetsize; ++i)
        for (j = 0; j < abetsize; ++j) {
            qij = fq1[i] * fq2[j] * exp(*Lambda * mtx[i][j]);
            sum += qij;
            *exp_res_score += qij * mtx[i][j];
        }

    exp_MSP_score = (log((double)n * n) + log(K)) / *Lambda;
    retval = (int)(exp_MSP_score / *exp_res_score + 0.5);

    *out_K = K;
    *out_H = H;
    *out_exp_msp = exp_MSP_score;

    free(prob); free(n1); free(n2); free(fq1); free(fq2);
    return retval;
}

/* ---- ASCII-to-binary tables, from dotplot.c:55 and dotplot.c:103 ---- */

static int atob_0[256] = {
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,23,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
    14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,
    NR, 0,20, 4, 3, 6,13, 7, 8, 9,NR,11,10,12, 2,NR,
    14, 5, 1,15,16,NR,19,17,22,18,21,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,
    NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR,NR
};

#define NN 5
static int ntob[256] = {
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
    NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN, 0,NN, 1,NN,NN,NN, 2,NN,NN,NN,NN,NN,NN, 4,NN,
    NN,NN,NN,NN, 3,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,
    NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN,NN
};

/* ---- BLOSUM62, verbatim from dotter.c:354 ---- */
static gint32 BLOSUM62[24][24] = {
    {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4 },
    { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4 },
    { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4 },
    { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
    {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4 },
    { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4 },
    { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
    {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4 },
    { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4 },
    { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4 },
    { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4 },
    { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4 },
    { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4 },
    { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4 },
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4 },
    {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4 },
    {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4 },
    { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4 },
    { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4 },
    {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4 },
    { -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4 },
    { -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4 },
    {  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4 },
    { -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1 }
};

/* ---- DNA matrix from dotter.c:2800; only the 4x4 corner matters ---- */
static gint32 DNAmtx[24][24];

static void make_dna_matrix(void) {
    int i, j;
    for (i = 0; i < 24; i++) for (j = 0; j < 24; j++) DNAmtx[i][j] = 0;
    for (i = 0; i < 6; i++)
        for (j = 0; j < 6; j++)
            DNAmtx[i][j] = (i < 4 && j < 4) ? (i == j ? 5 : -4) : -4;
}

/* ---- Fixture runner ---- */

static void run_fixture(const char *name, gint32 (*mtx)[24], int *tob,
                        int abetsize, const char *q, const char *s) {
    double lambda = 0, exp_res = 0, K = 0, H = 0, exp_msp = 0;
    int win = winsizeFromlambdak((gint32 (*)[24])mtx, tob, abetsize, q, s,
                                 &exp_res, &lambda, &K, &H, &exp_msp);
    /* Apply C-dotter's clamp [3, 50] in winsizeFromlambdak's caller. */
    int predicted = win;
    int window = win < 3 ? 3 : (win > 50 ? 50 : win);
    printf("%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\n",
           name, lambda, K, H, exp_res, exp_msp, predicted, window);
}

int main(void) {
    make_dna_matrix();

    /* Headers — column legend is documented at top of file. */
    printf("# name\tlambda\tK\tH\texp_res\texp_msp\tpredicted\twindow\n");

    /* === DNA fixtures (abetsize = 4) === */
    const char *dna_q1 = "ACGTACGTACGTACGTACGT";
    const char *dna_s1 = "GTGTACGAGCATCGTCTACT";
    /* The Rust port repeats the seed 40x to match its unit-test inputs. */
    static char dna_q1_long[20*40 + 1], dna_s1_long[20*40 + 1];
    for (int i = 0; i < 40; i++) {
        memcpy(dna_q1_long + i*20, dna_q1, 20);
        memcpy(dna_s1_long + i*20, dna_s1, 20);
    }
    dna_q1_long[20*40] = 0; dna_s1_long[20*40] = 0;
    run_fixture("dna_uniform", DNAmtx, ntob, 4, dna_q1_long, dna_s1_long);

    const char *dna_q2 = "AAAACCCGGGTTTAACAGCTAGCTACGATCGATCGATCGTAGCTAGCTAGCT";
    const char *dna_s2 = "TTTTGGGCCCAATTGCTAGCTACGATCGATCGATCGTAGCTAGCTAGCTACG";
    run_fixture("dna_at_rich_vs_gc", DNAmtx, ntob, 4, dna_q2, dna_s2);

    const char *dna_q3 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    /* Self-comparison */
    run_fixture("dna_self_repeat", DNAmtx, ntob, 4, dna_q3, dna_q3);

    /* === Protein fixtures (abetsize = 24) === */
    const char *prot_q1 =
        "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL"
        "RSALEFNAHVDEMVRLRREVGNQLEELQNRLREYIQRDHRGHEALQQYRVKQVHLDQEEIA";
    const char *prot_s1 =
        "MAATKRIIRQRYTIKHYVTRLREHIDHEEQVRKDLDEHKHRADRMLEELAGAILAAEHRLRD"
        "AREAFEQLLDKLEEHLRYAEELQEKFAKLERELAEHRLEEIEGRLAQAEEEFVEQHRRLENEL";
    run_fixture("prot_uniform", BLOSUM62, atob_0, 24, prot_q1, prot_s1);

    /* Mostly hydrophobic */
    const char *prot_q2 =
        "MFLLAAVAILALAIVAFLGLAILAVFAGLAAVILSFAGLLAVIFVASLLAGAILAVAFGLLI";
    const char *prot_s2 =
        "AILAFFGLAVILFAVALILGAVILAFAGAILVAFLLAVGAILSFAGLAILVAFAVALILSFG";
    run_fixture("prot_hydrophobic", BLOSUM62, atob_0, 24, prot_q2, prot_s2);

    /* Mostly charged */
    const char *prot_q3 =
        "MDDEERKKDEEDKKRRDEEKKRDEEEDKKRRDEEKKRDEEEDKKRRDEEKKRDEEEDKKRRDE";
    const char *prot_s3 =
        "EEDDKKRRDDEEKKRDDEEKKRRDDEEKKRRDDEEKKRDDEEKKRRDDEEKKRRDDEEKKRDD";
    run_fixture("prot_charged", BLOSUM62, atob_0, 24, prot_q3, prot_s3);

    return 0;
}
