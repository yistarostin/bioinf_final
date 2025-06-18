/*
Z-HUNT-2 computer program, Sept.19, 1990
*/

/*
serialized i/o to allow to run against large datasets
Updated i/o permissions to match posix. campt 1/10/2000
*/

/*
Turbo C compiler Ver. 1.5
Compact model, 8086/80286 instruction set, math emulator/8087/80287
Speed option, unsigned char
Run on IBM PC XT/AT or compatibles
*/

/*
Written by Ping-jung Chou, under the instruction of Pui S. Ho, according to the
paper "A computer aided thermodynamic approach for predicting the formation of
Z-DNA in naturally occurring sequences",
The EMBO Journal, Vol.5, No.10, pp2737-2744, 1986
With 0.22 kcal/mol/dinuc for mCG (Zacharias et al, Biochemistry, 1988, 2970)
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

/* creates a temporary file and mmaps the sequence into it */
#ifdef USE_MMAP
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#ifndef MAX_WINDOW_SIZE
#define MAX_WINDOW_SIZE 12
#endif
#ifndef THREADS
#include <sys/sysinfo.h>
#define THREADS get_nprocs_conf()
#endif

double linear_search(double x1, double x2, double tole, double (*func)(double));
double delta_linking(double dl);
double delta_linking_slope(double dl);
void show_probability(unsigned seqlength, float *dl, float *slope, float *probability,
                      char *sequence, char **antisyn, char *filename, int low, int high);

double linear_search(double x1, double x2, double tole, double (*func)(double)) {
    double f, fmid, dx, xmid, x;

    f = func(x1);
    fmid = func(x2);
    if (f * fmid >= 0.0)
        return x2;
    x = (f < 0.0) ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
    do {
        dx *= 0.5;
        xmid = x + dx;
        fmid = func(xmid);
        if (fmid <= 0.0)
            x = xmid;
    } while (fabs(dx) > tole);
    return x;
}

thread_local int terms;
double bztwist[MAX_WINDOW_SIZE];
thread_local double logcoef[MAX_WINDOW_SIZE], exponent[MAX_WINDOW_SIZE];
double _k_rt = -0.2521201;  /* -1100/4363 */
double sigma = 16.94800353; /* 10/RT */
thread_local double deltatwist;
double explimit = -600.0;

double delta_linking(double dl) {
    double sump, sumq, z, expmini;
    int i;

    expmini = 0.0;
    for (i = 0; i < terms; i++) {
        z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini)
            expmini = z;
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
    sump = sumq = 0.0;
    for (i = 0; i < terms; i++) {
        z = exp(exponent[i] + expmini);
        sumq += z;
        sump += bztwist[i] * z;
    }
    sumq += exp(_k_rt * dl * dl + sigma + expmini);
    return deltatwist - sump / sumq;
}

double delta_linking_slope(double dl) {
    double sump, sump1, sumq, sumq1, x, y, z, expmini;
    int i;

    expmini = 0.0;
    for (i = 0; i < terms; i++) {
        z = dl - bztwist[i];
        exponent[i] = z = logcoef[i] + _k_rt * z * z;
        if (z < expmini)
            expmini = z;
    }
    expmini = (expmini < explimit) ? explimit - expmini : 0.0;
    sump = sump1 = sumq = sumq1 = 0.0;
    x = 2.0 * _k_rt;
    for (i = 0; i < terms; i++) {
        z = dl - bztwist[i];
        y = exp(exponent[i] + expmini);
        sumq += y;
        sump += bztwist[i] * y;
        y *= z * x;
        sumq1 += y;
        sump1 += bztwist[i] * y;
    }
    y = exp(_k_rt * dl * dl + sigma + expmini);
    sumq += y;
    sumq1 += x * dl * y;
    return (sump1 - sump * sumq1 / sumq) / sumq;
} /* slope at delta linking = dl */

/* Delta BZ Energy of Dinucleotide */
double dbzed[4][17] = {
    /* AS-AS */
    {4.40, 6.20, 3.40, 5.20, 2.50, 4.40, 1.40, 3.30, 3.30, 5.20, 2.40, 4.20, 1.40, 3.40, 0.66, 2.40,
     4.26},
    /* SA-SA */
    {4.40, 2.50, 3.30, 1.40, 6.20, 4.40, 5.20, 3.40, 3.40, 1.40, 2.40, 0.66, 5.20, 3.30, 4.20, 2.40,
     4.26},
    /* AS-SA */
    {6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00,
     4.26},
    /* SA-AS */
    {6.20, 6.20, 5.20, 5.20, 6.20, 6.20, 5.20, 5.20, 5.20, 5.20, 4.00, 4.00, 5.20, 5.20, 4.00, 4.00,
     4.26}};

double expdbzed[4][17];                    /* exp(-dbzed/rt) */
thread_local int bzindex[MAX_WINDOW_SIZE]; /* dinucleotides */

int user_regret(void);
FILE *open_file(int mode, const char *filename, const char *typestr);
void assign_bzenergy_index(int nucleotides, char seq[]);
void best_anti_syn(int dinucleotides, float esum);
void anti_syn_energy(int dinculeotide, int dinucleotides, float esum);
unsigned input_sequence(FILE *file, int nucleotides, int showfile);
double assign_probability(double dl);

void assign_bzenergy_index(int nucleotides, char seq[]) {
    int i, j, idx;
    char c1, c2;

    i = j = 0;
    do {
        c1 = seq[i++];
        c2 = seq[i++];
        switch (c1) {
            case 'a':
                switch (c2) {
                    case 'n':
                        idx = 16;
                        break;
                    case 'a':
                        idx = 0;
                        break;
                    case 't':
                        idx = 1;
                        break;
                    case 'g':
                        idx = 2;
                        break;
                    case 'c':
                        idx = 3;
                        break;
                }
                break;
            case 't':
                switch (c2) {
                    case 'n':
                        idx = 16;
                        break;
                    case 'a':
                        idx = 4;
                        break;
                    case 't':
                        idx = 5;
                        break;
                    case 'g':
                        idx = 6;
                        break;
                    case 'c':
                        idx = 7;
                        break;
                }
                break;
            case 'g':
                switch (c2) {
                    case 'n':
                        idx = 16;
                        break;
                    case 'a':
                        idx = 8;
                        break;
                    case 't':
                        idx = 9;
                        break;
                    case 'g':
                        idx = 10;
                        break;
                    case 'c':
                        idx = 11;
                        break;
                }
                break;
            case 'c':
                switch (c2) {
                    case 'n':
                        idx = 16;
                        break;
                    case 'a':
                        idx = 12;
                        break;
                    case 't':
                        idx = 13;
                        break;
                    case 'g':
                        idx = 14;
                        break;
                    case 'c':
                        idx = 15;
                        break;
                }
                break;
            case 'n':
                switch (c2) {
                    case 'a':
                        idx = 16;
                        break;
                    case 't':
                        idx = 16;
                        break;
                    case 'g':
                        idx = 16;
                        break;
                    case 'n':
                        idx = 16;
                        break;
                    case 'c':
                        idx = 16;
                        break;
                }
                break;
        }
        bzindex[j++] = idx;
    } while (i < nucleotides);
}

thread_local double bzenergy[MAX_WINDOW_SIZE], best_bzenergy[MAX_WINDOW_SIZE]; /* dinucleotides */
thread_local float best_esum; /* assigned before call to anti_syn_energy() */
thread_local char best_antisyn[2 * MAX_WINDOW_SIZE + 1],
    antisyn[2 * MAX_WINDOW_SIZE + 1]; /* nucleotides */

void best_anti_syn(int dinucleotides, float esum) {
    int i;
    double dl, slope;

    if (esum < best_esum) {
        best_esum = esum;
        for (i = 0; i < dinucleotides; i++)
            best_bzenergy[i] = bzenergy[i];
        strcpy(best_antisyn, antisyn);
    }
}

void anti_syn_energy(int din, int dinucleotides, float esum) {
    int i, nucleotides;
    float e;

    nucleotides = 2 * din;

    antisyn[nucleotides] = 'A';
    antisyn[nucleotides + 1] = 'S';
    i = (din == 0) ? 0 : ((antisyn[nucleotides - 1] == 'S') ? 0 : 3);
    e = dbzed[i][bzindex[din]];
    esum += e;
    bzenergy[din] = expdbzed[i][bzindex[din]];
    if (++din == dinucleotides)
        best_anti_syn(dinucleotides, esum);
    else
        anti_syn_energy(din, dinucleotides, esum);
    esum -= e;
    din--;

    antisyn[nucleotides] = 'S';
    antisyn[nucleotides + 1] = 'A';
    i = (din == 0) ? 1 : ((antisyn[nucleotides - 1] == 'A') ? 1 : 2);
    esum += dbzed[i][bzindex[din]];
    bzenergy[din] = expdbzed[i][bzindex[din]];
    if (++din == dinucleotides)
        best_anti_syn(dinucleotides, esum);
    else
        anti_syn_energy(din, dinucleotides, esum);
}

thread_local char tempstr[128];
char *sequence;
#ifdef USE_MMAP
int sequencefile;
#endif

int user_regret(void) {
    char c;
    do {
        fgets(tempstr, 128, stdin);
        c = tempstr[0];
    } while (c == 0);
    return (c == '@') ? 1 : 0;
}

FILE *open_file(int mode, const char *filename, const char *typestr) {
    static const char *iostr[] = {"output", "input"};
    static const char *rwstr[] = {"w", "r"};
    char *fullfile;

    FILE *file;
    file = NULL;
    fullfile = (char *)malloc(sizeof(char) * (strlen(filename) + strlen(typestr) + 2));
    strcpy(fullfile, filename);
    if (strlen(typestr) != 0) {
        strcat(fullfile, ".");
        strcat(fullfile, typestr);
    }
    printf("opening %s\n", fullfile);

    file = fopen(fullfile, rwstr[mode]);
    free(fullfile);
    return file;
}

unsigned input_sequence(FILE *file, int nucleotides, int showfile) {
    unsigned length, i, j;
    char c;
#ifdef USE_MMAP
    FILE *OUTPUT;
#endif

    printf("inputting sequence\n");

    length = 0; /* count how many bases */
    while (j = 0, fgets(tempstr, 128, file) != NULL) {
        while ((c = tempstr[j++]) != 0) {
            if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n' || c == 'A' || c == 'T' ||
                c == 'G' || c == 'C' || c == 'N') {
                length++;
            }
        }
    }

#ifndef USE_MMAP
    sequence = (char *)malloc(length + nucleotides);
#else
    sequence = (char *)malloc(nucleotides);
    OUTPUT = (FILE *)fopen("/usr/local/apache/htdocs/zhunt/temp", "w");
#endif

    rewind(file);

    if (showfile) {
        printf("\n");
    }
    i = 0;
    while (j = 0, fgets(tempstr, 128, file) != NULL) {
        while ((c = tempstr[j++]) != 0) {
            if (c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'n' || c == 'A' || c == 'T' ||
                c == 'G' || c == 'C' || c == 'N') {
                if (c == 'N') {
                    c = 'n';
                }
                if (c == 'A') {
                    c = 'a';
                }
                if (c == 'T') {
                    c = 't';
                }
                if (c == 'G') {
                    c = 'g';
                };
                if (c == 'C') {
                    c = 'c';
                }
#ifndef USE_MMAP
                sequence[i++] = c;
#else
                fprintf(OUTPUT, "%c", c);
                if (i < nucleotides) {
                    sequence[i++] = c;
                }
#endif
            }
        }
    }
    for (j = 0; j < nucleotides; j++) /* assume circular nucleotides */
    {
#ifndef USE_MMAP
        sequence[i++] = sequence[j];
#else
        fprintf(OUTPUT, "%c", sequence[j]);
#endif
    }
#ifdef USE_MMAP
    close(OUTPUT);
    sequencefile = open("/usr/local/apache/htdocs/zhunt/temp", O_RDWR, NULL);
    free(sequence);
    sequence = (char *)mmap(0, length, PROT_READ | PROT_WRITE, MAP_SHARED, sequencefile, 0);
#endif

    return length;
}

/* calculate the probability of the value 'dl' in a Gaussian distribution */
/* from "Data Reduction and Error Analysis for the Physical Science" */
/* Philip R. Bevington, 1969, McGraw-Hill, Inc */

double assign_probability(double dl) {
    static double average = 29.6537135;
    static double stdv = 2.71997;
    static double _sqrt2 = 0.70710678118654752440; /* 1/sqrt(2) */
    static double _sqrtpi = 0.564189583546;        /* 1/sqrt(pi) */

    double x, y, z, k, sum;

    z = fabs(dl - average) / stdv;
    x = z * _sqrt2;
    y = _sqrtpi * exp(-x * x);
    z *= z;
    k = 1.0;
    sum = 0.0;
    do {
        sum += x;
        k += 2.0;
        x *= z / k;
    } while (sum + x > sum);
    z = 0.5 - y * sum; /* probability of each tail */
    return (dl > average) ? z : 1.0 / z;
}

double find_delta_linking(int dinucleotides);
void calculate_zscore(double a, int maxdinucleotides, int min, int max, char *filename);

void print_array(int points[], unsigned from, unsigned to, int printwidth);
void soft_copy(int width, int widthbyte, char **screen);

void generate_random_sequence(int dinucleotides, char seq[]);
void initiate_random(long seed);
double uniform_random(void);
void random_distribution(double a, int maxdinucleotides, char *filename);
void run_distribution(double a, int maxdinucleotides, FILE *disfile);
void print_distribution(char *disname, FILE *disfile);

double find_delta_linking(int dinucleotides) {
    double sum;
    int i, j;

    for (i = 0; i < dinucleotides; i++)
        bzenergy[i] = 1.0;
    for (i = 0; i < dinucleotides; i++) {
        sum = 0.0;
        for (j = 0; j < dinucleotides - i; j++) {
            bzenergy[j] *= best_bzenergy[i + j];
            sum += bzenergy[j];
        }
        logcoef[i] = log(sum);
    }
    terms = dinucleotides;
    return linear_search(10.0, 50.0, 0.001, delta_linking);
}

static void calc(unsigned i, int nucleotides, unsigned fromdin, unsigned todin, float initesum,
                 double a, FILE *file, pthread_mutex_t *output_lock) {
    static double pideg = 57.29577951; /* 180/pi */

    assign_bzenergy_index(nucleotides, sequence + i);
    double bestdl = 50.0;
    char bestantisyn[nucleotides + 1];
    for (unsigned din = fromdin; din <= todin; din++) {
        best_esum = initesum;
        deltatwist = a * (double)din;
        antisyn[2 * din] = 0;
        anti_syn_energy(0, din, 0.0); /* esum = 0.0 */
        double dl = find_delta_linking(din);
        if (dl < bestdl) {
            bestdl = dl;
            strcpy(bestantisyn, best_antisyn);
        }
    }
    unsigned j = strlen(bestantisyn);
#ifndef PROB_ONLY
    double slope = atan(delta_linking_slope(bestdl)) * pideg;
#endif
    double probability = assign_probability(bestdl);
    if (probability < 300)
        return;

    pthread_mutex_lock(output_lock);

#ifndef PROB_ONLY
    fprintf(file, "%i %i %i %7.3lf %7.3lf %le ", i + 1, i + 1 + j, j, bestdl, slope, probability);
#else
    fprintf(file, " %7.3lf %le\n", bestdl, probability);
#endif
    for (unsigned k = 0; k < j; k++)
        fprintf(file, "%c", sequence[i + k]);
    fprintf(file, "   %s\n", bestantisyn);
    fflush(file);

    pthread_mutex_unlock(output_lock);
}

struct Params {
    int threads;
    unsigned seqlength;
    int nucleotides;
    unsigned fromdin;
    unsigned todin;
    float initesum;
    double a;
    FILE *file;
} static params;

struct ThreadContext {
    size_t thread_id;
    pthread_mutex_t *output_lock;
};

static int run_thread(void *ctx) {
    struct ThreadContext *thread_ctx = (struct ThreadContext *)ctx;
    for (size_t i = thread_ctx->thread_id; i < params.seqlength; i += params.threads) {
        calc(i, params.nucleotides, params.fromdin, params.todin, params.initesum, params.a,
             params.file, thread_ctx->output_lock);
    }
    return 0;
}

void calculate_zscore(double a, int maxdinucleotides, int min, int max, char *filename) {
    static double pideg = 57.29577951; /* 180/pi */
    char *bestantisyn;
    FILE *file;
    unsigned seqlength, i, j, k;
    int fromdin, todin, din, nucleotides;
    long begintime, endtime;
    double dl, slope, probability, bestdl;
    float initesum;

    fromdin = min;
    todin = max;
    printf("calculating zscore\n");

    file = open_file(1, filename, "");
    if (file == NULL) {
        printf("couldn't open %s!\n", filename);
        return;
    }
    seqlength = input_sequence(file, 2 * maxdinucleotides, 0);
    fclose(file);

    file = open_file(0, filename, "Z-SCORE");
    if (file == NULL) {
#ifndef USE_MMAP
        free(sequence);
#else
        munmap(sequence, seqlength);
        close(sequencefile);
#endif
        return;
    }

    if (todin > maxdinucleotides) {
        todin = maxdinucleotides;
    }
    if (fromdin > todin) {
        fromdin = todin;
    }
    nucleotides = 2 * todin;

#ifndef PROB_ONLY
    fprintf(file, "%s %u %d %d\n", filename, seqlength, fromdin, todin);
#endif

    a /= 2.0;
    initesum = 10.0 * todin;
    bestantisyn = (char *)malloc(nucleotides + 1);

    time(&begintime);

    params.threads = THREADS;
    params.seqlength = seqlength;
    params.nucleotides = nucleotides;
    params.fromdin = fromdin;
    params.todin = todin;
    params.initesum = initesum;
    params.a = a;
    params.file = file;

    pthread_mutex_t output_lock;
    pthread_mutex_init(&output_lock, NULL);

    struct ThreadContext ctxs[params.threads];

    thrd_t threads[params.threads];
    for (size_t t = 0; t < params.threads; ++t) {
        ctxs[t] = (struct ThreadContext){
            .thread_id = t,
            .output_lock = &output_lock,
        };
        thrd_create(&threads[t], run_thread, (void *)&ctxs[t]);
    }
    for (size_t t = 0; t < params.threads; ++t) {
        thrd_join(threads[t], NULL);
    }

    pthread_mutex_destroy(&output_lock);

    time(&endtime);

    free(bestantisyn);
    fclose(file);
    printf("\n run time=%ld sec\n", endtime - begintime);
#ifndef USE_MMAP
    free(sequence);
#else
    munmap(sequence, seqlength);
    close(sequencefile);
#endif
}

void analyze_zscore(char *filename, int low, int high) {
    float *dl, *slope, *probability;
    unsigned seqlength, i, select;
    char **antisyn;
    FILE *file;
    int fromdin, todin, nucleotides;

    printf("use min/max %d %d\n", low, high);
    printf("analyzing_zscore\n");

    file = open_file(1, filename, "Z-SCORE");
    if (file == NULL) {
        printf("couldn't open %s.Z-SCORE!\n", filename);
        return;
    }
    fscanf(file, "%s %u %d %d", tempstr, &seqlength, &fromdin, &todin);
    nucleotides = 2 * todin;
    dl = (float *)calloc(seqlength, sizeof(float));
    slope = (float *)calloc(seqlength, sizeof(float));
    probability = (float *)calloc(seqlength, sizeof(float));
    antisyn = (char **)calloc(seqlength, sizeof(char *));

    for (i = 0; i < seqlength; i++) {
        fscanf(file, "%f %f %f %s", dl + i, slope + i, probability + i, tempstr);
        antisyn[i] = strdup(tempstr);
    }
    fclose(file);

    file = open_file(1, filename, "");
    input_sequence(file, nucleotides, 0);
    fclose(file);
    /* uncomment to print probability file
     show_probability( seqlength, dl, slope, probability, sequence, antisyn,filename,low,high);
     */
    free(dl);
    free(slope);
    free(probability);
    for (i = 0; i < seqlength; i++)
        free(antisyn[i]);
    free(antisyn);

#ifndef USE_MMAP
    free(sequence);
#else
    munmap(sequence, seqlength);
    close(sequencefile);
#endif
}

long RAm = 714025, RAa = 1366, RAc = 150889, RAseed = 0, RAy, RAr[97];

void initiate_random(long seed) {
    int i;

    seed = (RAc - seed) % RAm;
    if (seed < 0)
        seed = -seed;
    for (i = 0; i < 97; i++)
        RAr[i] = seed = (RAa * seed + RAc) % RAm;
    seed = (RAa * seed + RAc) % RAm;
    RAseed = RAy = seed;
}

double uniform_random(void) {
    int i;

    i = (double)97.0 * RAy / (double)RAm;
    RAy = RAr[i];
    RAr[i] = RAseed = (RAa * RAseed + RAc) % RAm;
    return (double)RAy / RAm;
}

void generate_random_sequence(int dinucleotides, char seq[]) {
    static char rs1[] = "aaaattttggggcccc";
    static char rs2[] = "atgcatgcatgcatgc";
    int i, j, r;

    for (i = j = 0; i < dinucleotides; i++) {
        r = uniform_random() * 17.0;
        seq[j++] = rs1[r];
        seq[j++] = rs2[r];
    }
}

void random_distribution(double a, int maxdinucleotides, char *filename) {
    int select;
    FILE *file;

    do {
        printf("\n 0. Stop");
        printf("\n 1. Run distribution");
        printf("\n\n Select: ");
        scanf("%d", &select);
        if (select > 0) {
            file = open_file(select - 1, filename, "DISTRIBUTION");
            if (file == NULL)
                continue;
            if (select == 1)
                run_distribution(a, maxdinucleotides, file);
            fclose(file);
        }
    } while (select != 0);
}

void run_distribution(double a, int maxdinucleotides, FILE *disfile) {
    int i, nucleotides, dinucleotides, levels;
    long begintime, endtime, repeat, k, *distribution;
    float initesum;
    double dl, sumdl, sumdldl, dlaverage, dlstdv, dlfrom, dlto, dlstep;

    printf("\n Window size (dinucleotides) and sample size: ");
    if (user_regret())
        return;
    sscanf(tempstr, "%d %ld", &dinucleotides, &repeat);
    if (dinucleotides > maxdinucleotides)
        dinucleotides = maxdinucleotides;
    nucleotides = 2 * dinucleotides;

    printf("\n FROM, TO & STEP linking differences for distribution: ");
    if (user_regret())
        return;
    sscanf(tempstr, "%lf %lf %lf", &dlfrom, &dlto, &dlstep);

    sequence = (char *)malloc(nucleotides);
    /* half of the window convert from B to Z */
    deltatwist = a * (double)dinucleotides / 2.0;
    initesum = 10.0 * dinucleotides; /* for best_anti_syn() use */

    levels = (dlto - dlfrom) / dlstep + 1;
    distribution = (long *)calloc(levels, sizeof(long));
    for (i = 0; i < levels; i++)
        distribution[i] = 0;
    sumdl = sumdldl = 0.0;

    time(&begintime);
    initiate_random(begintime);

    for (k = 0; k < repeat; k++) {
        generate_random_sequence(dinucleotides, sequence);
        assign_bzenergy_index(nucleotides, sequence);
        best_esum = initesum;
        anti_syn_energy(0, dinucleotides, 0.0);
        dl = find_delta_linking(dinucleotides);
        sumdl += dl;
        sumdldl += dl * dl;
        i = (dl - dlfrom) / dlstep;
        if (i < 0)
            i = 0;
        if (i > levels - 1)
            i = levels - 1;
        distribution[i]++;
    }
    time(&endtime);
    dlaverage = sumdl / repeat;
    dlstdv = sqrt((sumdldl - sumdl * dlaverage) / repeat);

    fprintf(disfile, " %d %d %ld %lf %lf %lg %lg %lg\n", dinucleotides, levels, repeat, dlaverage,
            dlstdv, dlfrom, dlto, dlstep);
    for (i = 0; i < levels; i++)
        fprintf(disfile, "%ld\n", distribution[i]);

    printf("\n run time: %ld sec\n", endtime - begintime);
    free(distribution);
    free(sequence);
}

void show_probability(unsigned seqlength, float *dl, float *slope, float *probability,
                      char *sequence, char **antisyn, char *filename, int low, int high) {
    unsigned i, j, k;
    FILE *file;

    printf("show_probability\n");
    printf("min/max %d %d\n", low, high);

    file = open_file(0, filename, "probability");

    if (file == NULL) {
        printf("couldn't open %s.probability!\n", filename);
        return;
    }
    fprintf(file, "%s %u %i %i\n", filename, seqlength, low, high);
    for (i = 0; i < seqlength; i++) {
        j = strlen(antisyn[i]);
        fprintf(file, " %i %i %i %7.3f %7.3f  %10.3e  ", i + 1, i + 1 + j, j, dl[i], slope[i],
                probability[i]);

        for (k = 0; k < j; k++)
            fprintf(file, "%c", sequence[i + k]);
        fprintf(file, "   %s\n", antisyn[i]);
    }

    fclose(file);
}

int main(int argc, char *argv[]) {
    static double rt = 0.59004;       /* 0.00198*298 */
    static double a = 0.357, b = 0.4; /* a = 2 * (1/10.5 + 1/12) */
    double ab;
    int i, j, nucleotides, dinucleotides, select;
    int min, max;
    int low, high;

    if (argc < 5) {
        printf("usage: zhunt windowsize minsize maxsize datafile\n");
        exit(1);
    }
    dinucleotides = atoi((char *)argv[1]);
    min = atoi((char *)argv[2]);
    max = atoi((char *)argv[3]);
    low = atoi((char *)argv[2]);
    high = atoi((char *)argv[3]);

    if (dinucleotides > MAX_WINDOW_SIZE) {
        fprintf(stderr, "windows size is greater than maximum: %d > %d\n", dinucleotides,
                MAX_WINDOW_SIZE);
        return 1;
    }

    printf("dinucleotides %d\n", dinucleotides);
    printf("min/max %d %d\n", min, max);
    printf("min/max %d %d\n", low, high);
    printf("operating on %s\n", (char *)argv[4]);

    nucleotides = 2 * dinucleotides;

    ab = b + b;
    for (i = 0; i < dinucleotides; i++) {
        ab += a;
        bztwist[i] = ab;
    }

    for (i = 0; i < 4; i++)
        for (j = 0; j < 17; j++)
            expdbzed[i][j] = exp(-dbzed[i][j] / rt);

    calculate_zscore(a, dinucleotides, min, max, (char *)argv[4]);
#ifndef PROB_ONLY
    analyze_zscore((char *)argv[4], low, high);
#endif

    return 0;
}
