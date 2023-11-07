/* Cluster identification routine based on ten Wolde's local orientational
 * order parameter analysis
 */





//    LIBRARIES
// ===============
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <math.h>               // need to compile with `-lm` flag
#include <string.h>
#include <getopt.h>





//    PREPROCESSOR CONSTANTS
// ============================
#define INVFPI 0.07957747154594766788444188168625718101 // = 1 / (4 * pi)
#define MAX_PART_CELL 40
#define MAX_NEIGHBORS 50





//    MACROS
// ============
#define OPTIONAL_ARGUMENT_IS_PRESENT \
    ((optind < argc-1 && optarg == NULL && optind < argc && argv[optind][0] != '-') \
     ? (uintptr_t) (optarg = argv[optind++]) \
     : (optarg != NULL))





//    PROTOTYPES
// ================
// structures
// ----------------
typedef struct Particle Particle;
typedef struct compl compl;
typedef struct bndT bndT;
typedef struct blistT blistT;
typedef struct posnb posnb;
typedef struct clusterBook clusterBook;
// functions
// ----------------
int parse_input(int argc, char* argv[]);
void initReading(FILE* initfile);
void readCoords(FILE* initfile);
void writeCoords(char *filename, int fluidlike, int snap);
int buildCL(int usecase);
int retrieveIndex(int u, int v, int w);
int findClusters(void);
void build_nblist(int method);
int compare(const void * a, const void * b);
compl* calc_order(void);
double dotprod(compl *vec1, compl *vec2, int l);
double sqr(double x);
void compute_order(int l, bndT *bnd, compl *res1, compl *res2);
float plgndr(int l, int m, double x);
double facs(int l, int m);
double gammln(double xx);
double minpow(int m);
int* calc_conn(compl* orderp);
void setcluss(int pn, compl* orderp, int cn);
int calc_clusters(int* conn, compl* orderp);
void saveLogBook(void);
void saveDGndata(void);





//    GLOBAL VARIABLES
// ======================
// hard spheres system
// ----------------------
int N = 2000;                           // [READ FROM IMPUT] number of particles in the system
double boxx = 0.0f;
double boxy = 0.0f;
double boxz = 0.0f;
double sigma = 1.0f;                    // unit of length = 1 = max particle size
Particle *particles = NULL;             // [READ FROM INPUT] pointer to the table of particles data
int PBC = 1;                            // [TUNABLE] [-p] choice to use periodic boundary conditions for analysis, default is yes (1)
// cell lists
// -----------------------
Particle **CLTable = NULL;
double sCellx = 0.0f;
int nCellx = 0;
double sCelly = 0.0f;
int nCelly = 0;
double sCellz = 0.0f;
int nCellz = 0;
int nCell = 0;
// bop
// ----------------------
int* conn = NULL;
int* connections = NULL;
compl* order = NULL;
compl* orderp = NULL;
int* cluss = NULL;
int* size = NULL;
double bndLength = 1.4f;                // [TUNABLE] [-r] distance cutoff for bonds, if used
double bnd_cutoff = 0.7f;               // [TUNABLE] [-b] dot product cutoff to be call a connected particle
                                        // `c` in ten Wolde's work (`c=0.5`), `d_c` in Filion's work (`d_c=0.7`)
int nbnd_cutoff = 4;                    // [TUNABLE] [-n] minimum number of connected particles to be called a crystaline particle
                                        // `nc` in ten Wolde's work (`nc=?`), `xi_c` in Filion's work (`4<xi_c<10`)
double obnd_cutoff = 0.0f;              // [TUNABLE] [-o] dot product cutoff to be in the same cluster
                                        // additional criterion for cluster identification: =0 for all touching clusters, =0.9 to see defects
double maxr2 = 0.0f;
blistT *blist;
int NNidMethod = 0;                     // [TUNABLE] choice of method for NN identification
clusterBook* logBook = NULL;            // bookkeeping table for the clusters identified in the system
int cs;                                 //cluster size
// ensemble average
// ----------------------
double* Nn = NULL;                      // table for calculating the free energy of formation of a cluster of size n
// input/output
// ----------------------
char* input_filename;                   // name of the initial configuration input file
char movie_filename[200];               // name of the movie output file
char clusterlog_filename[200];          // name of the clusters log output file
char DGndata_filename[200];             // name of the DGn data output file
long cursor_end;                        // position of the EOF for the input SOURCE file
long cursor_current;                    // current reading position in the input SOURCE file
int input_filetype = 0;                 // [TUNABLE] [-f] input file type: 0 is for standard .sph, 1 is for .xyz datafiles
int output_filetype = 0;                // [TUNABLE] [-F] output file type: 0 is for standard .sph, 1 is for .xyz
int highlightAll = 1;                   // [TUNABLE] [-L] choice to highlight all found clusters (1) or largest found cluster only (0)
int snapshot_count = 0;                 // count for the number of read snapshots
int savelogs = 0;                       // [TUNABLE] [-l] choice to save logs
int energy = 0;                         // [TUNABLE] [-G] choice to calculate free energy to create a nucleus of size n
int save_movie = 1;                     // [TUNABLE] [-M] choice to save the movie of clusters, default is yes (1)
int snap_count = 0;                     // snapshot counter




//    STRUCTURES
// ================
struct Particle{
/* Structure:  Particle
 * --------------------
 * Particle in the system, implemented for the use of cell lists (CL)
 */
	double x;               // reduced x coordinate 
	double y;	            // reduced y coordinate
    double z;               // reduced z coordinate
	double r;	            // reduced radius --- redundant with sigma ... might remove later on
    char type;              // particle type --- for visualization purposes only
    Particle *next;         // pointer to the next particle in CL
    Particle *prev;         // pointer to the previous particle in CL
    int CLindex;            // index of CL in which particle is
    int index;              // index of the particle, in [0;N-1]
    int clusterindex;       // index of the cluster the particle belongs to; 0 if fluidlike
};



struct compl{
/* Structure:  compl
 * -----------------
 * Complex number
 */
    double re;              // real part
    double im;              // imaginary part
};



struct bndT{
/* Structure: bndT
 * ---------------
 * Geometrical information about the bond between two neighbouring particles
 */
    double nz;
    double si;
    double co;
    int n;
};



struct blistT{
/* Structure:  blistT
 * ------------------
 * List of particles connected to particles i
 */
    int n;                  // number of particles connected to particle i
    bndT * bnd;             // description of the bond connecting particle j to particle i
};



struct posnb{
/* Structure:   posnb
 * ------------------
 * List and info about the NN of a given partice
 */
    Particle* part;         // pointer to a NN particle
    double dist;            // distance between the two particles
    double dx;              // x-axis distance between the two particles
    double dy;              // y-axis distance between the two particles              
    double dz;              // z-axis distance between the two particles
};



struct clusterBook{
/* Structure:   clusterBook
 * ------------------------
 * Bookkeeping for the clusters identified in the system
 */
    int bookmark;           // number of cluster sizes identified
    int* size;              // size of the identified cluster
    int* qt;                // number of identified clusters of a given size
};





//    FUNCTIONS
// ===============
int parse_input(int argc, char* argv[]){
/* Function:    parse_input
 * ------------------------
 * Parser for the programme execution
 *
 * argc:                ?
 * argv:                ?
 *
 * return:      0 for normal termination
 */
        // Parsing of command line arguments
	struct option longopt[] = {
        {"bnd-cutoff",required_argument,NULL,'b'},
        {"nbnd-cutoff",required_argument,NULL,'n'},
        {"obnd-cutoff",required_argument,NULL,'o'},
        {"rc",required_argument,NULL,'r'},
        {"highlight-all",no_argument,NULL,'L'},
        {"input-filetype",required_argument,NULL,'f'},
        {"output-filetype",required_argument,NULL,'F'},
        {"NN-method",required_argument,NULL,'m'},
        {"save-logs",no_argument,NULL,'l'},
        {"DGn",no_argument,NULL,'G'},
        {"save-movie",required_argument,NULL,'M'},
        {"PBC",required_argument,NULL,'p'},
		{"help",no_argument,NULL,'h'},
		{0,0,0,0}
        };

        int c;
        while ((c = getopt_long(argc, argv, "b:n:o:r:Lf:F:m:lGM:p:h", longopt, NULL)) != -1){
                switch (c){
                        case 'b':
                                if (sscanf(optarg, "%lf", &bnd_cutoff) != 1){
                                        printf("[clusterID] Could not parse bnd_cutoff value.\n");
                                        exit(3);
                                }
                                break;
                        case 'n':
                                if (sscanf(optarg, "%d", &nbnd_cutoff) != 1){
                                        printf("[clusterID] Could not parse nbnd_cutoff value.\n");
                                        exit(3);
                                }
                                break;
                        case 'o':
                                if (sscanf(optarg, "%lf", &obnd_cutoff) != 1){
                                        printf("[clusterID] Could not parse obnd_cutoff value.\n");
                                        exit(3);
                                }
                                break;
                        case 'r':
                                if (sscanf(optarg, "%lf", &bndLength) != 1){
                                        printf("[clusterID] Could not parse bndLength value.\n");
                                        exit(3);
                                }
                                NNidMethod = 2;
                                break;
                        case 'L':
                                highlightAll = 0;
                                break;
                        case 'f':
                                if (sscanf(optarg, "%d", &input_filetype) != 1){
                                        printf("[clusterID] Could not parse input_filetype.\n");
                                        exit(3);
                                }
                                break;
                        case 'F':
                                if (sscanf(optarg, "%d", &output_filetype) != 1){
                                        printf("[clusterID] Could not parse output_filetype.\n");
                                        exit(3);
                                }
                                break;
                        case 'm':
                                if (sscanf(optarg, "%d", &NNidMethod) != 1){
                                        printf("[clusterID] Could not parse NNidMethod.\n");
                                        exit(3);
                                }
                                break;
                        case 'l':
                                savelogs = 1;
                                break;
                        case 'G':
                                energy = 1;
                                break;
                        case 'M':
                                if (sscanf(optarg, "%d", &save_movie) != 1) {
                                        printf("[clusterID] Could not parse choice to save cluster movie.\n");
                                        exit(3);
                                }
                        case 'p':
                                if (sscanf(optarg, "%d", &PBC) != 1) {
                                        printf("[clusterID] Could not parse choice to use PBC in analysis.\n");
                                        exit(3);
                                }
                                break;
                        case 'h':
                                printf("[clusterID]\n * Usage: ./clusterID [OPTIONS] SOURCE\n * Description: identifies clusters based on the choice of parameters / method and writes a .sph file were they are highlighted.\n * Options:\n * -b [double]: ten Wolde's local orientational order parameter dot-product cutoff for two particles to be identified as connected\n * -n [int]: number of connections cutoff for a particle to be identified as solid-like\n * -o [double]: ten Wolde's local orientational order parameter dot-product cutoff for two particles to be identified as part of the same cluster\n * -r [double]: radius cutoff for Nearest Neighbours (NN) identification. Automatically triggers choice of cutoff radius method for NN identification\n * -L: only highlights the largest identified cluster\n * -f [0 or 1]: choice of SOURCE file type among:\n *           0 for .sph files (default)\n *           1 for .xyz files\n * -F [0 or 1]: choice of OUTPUT file type, choice is the same as -f (default if 0 for .sph format)\n * -m [0 or 1 or 2]: choice of NN identification method among:\n *           0 for SANN method (default)\n *           1 for 12NN method\n *           2 for cutoff radius method\n * -l: saves cluster logs to a file\n * -G: calculates free energy for the formation of a cluster of size n, saves it to a file\n * -M [0 or 1]: choice to save the movie of clusters (default is 1)\n * -p [0 or 1]: use PBC for analysis (default is 1)\n * -h: displays this message and exit\n");
                                exit(0);
                }
        }

        // Wrong usage
        if(optind>argc-1){
		printf("[clusterID] Usage: ./clusterID [OPTIONS] SOURCE\n [clusterID] Help: ./clusterID -h\n");
		
        if (input_filetype == 1) PBC = 0;
        exit(0);
	}
        
        // Reading input SOURCE filename
        input_filename = argv[optind];
        
        return 0;
}



void initReading(FILE* initfile){
/* Function:    initReading
 * ------------------------
 * Initializes the reading of the supplied SOURCE file
 *
 * *initfile:   pointer to the supplied SOURCE file
 */
        // Localization of the EOF
        fseek(initfile, 0, SEEK_END);
        cursor_end = ftell(initfile);
        rewind(initfile);
        cursor_current = ftell(initfile);
}



int mygetline(char* str, FILE* f) {
        int comment = 1;
        while (comment) {
                if (!fgets(str, 255, f)) return -1;
                if (str[0] != '#') comment = 0;
        }
        return 0;
}


void readCoords(FILE* initfile){
/* Function:    readCoords
 * -----------------------
 * Initializes the table of data of all N particles by reading
 * it from the supplied SOURCE file; particles coordinates
 * are written in reduced units; also retrieves **cubic** box size
 *
 * *initfile:   pointer to the supplied SOURCE file
 */
        char buffer[255];
        int ftmp;
        // Initialization of the number of particles N
        switch (input_filetype){
                case 0:
                        mygetline(buffer,initfile);
                        ftmp = sscanf(buffer, "%d\n", &N);
                        if (ftmp != 1) {printf("Error!\n"); exit(3);}
                        break;
                case 1:
                        if (fscanf(initfile, "%d%*c", &N) != 1) {
                                printf("[findClusters] Error: could not parse input file. Check input filetype.\n");
                                exit(0);
                        }
                        break;
                default:
                        printf("[findClusters] Error: wrong input filetype specified.\n");
                        exit(0);
        }
        // Memory allocations
        particles = malloc(N * sizeof(*particles));
        conn = malloc(N * sizeof(*conn));
        connections = malloc(N * sizeof(*connections));
        order = malloc(N * (2 * 6 + 1) * sizeof(*order));
        orderp = malloc(N * (2 * 6 + 1) * sizeof(*orderp));
        cluss = malloc(N * sizeof(*cluss));
        size = malloc(N * sizeof(*size));
        blist = malloc(N * sizeof(*blist));
        for (int i = 0; i < N; i++)
                blist[i].bnd = malloc(MAX_NEIGHBORS * sizeof(bndT));
        logBook = malloc(sizeof(*logBook));
        logBook->size = malloc(N * sizeof(*logBook->size));
        logBook->qt = malloc(N * sizeof(*logBook->qt));
        // More memory allocation
        static int allocate = 1;
        if (allocate){
                Nn = malloc(N * sizeof(*Nn));
                memset(Nn, (double) 0.0f, N*sizeof(*Nn));
                allocate = 0;
        }

        // Read 2nd line
        switch (input_filetype){
                case 0:
                        if (fscanf(initfile, "%lf %lf %lf%*c", &boxx, &boxy, &boxz) != 3) {printf("error here 0\n"); exit(0);}
                        break;
                case 1:
                        if (fscanf(initfile, "%*s") != 0) {printf("error here 1\n"); exit(0);}
                        break;
        }

        // Populate table of particles data in standard units
        switch (input_filetype){
                case 0: ;
                        double rTemp = 0.0f;
                        for (int i = 0; i < N; i++){
                                if (fscanf( initfile,
                                            "%c %lf %lf %lf %lf%*c",
                                            &particles[i].type,
                                            &particles[i].x,
                                            &particles[i].y,
                                            &particles[i].z,
                                            &particles[i].r
                                          )
                                                != 5) {
                                        if (fscanf(initfile, "%c %lf  %lf  %lf  %lf%*c", &particles[i].type, &particles[i].x, &particles[i].y, &particles[i].z, &particles[i].r) != 5) {
                                                printf("error -1\n");
                                                exit(0);
                                        }
                                }
                                // Find larger radius to rescale everything to obtain
                                // the expected sigma = 1 for the (larger) particle size
                                if (particles[i].r > rTemp)
                                        rTemp = particles[i].r;
                        }
                        // Rescale everything w.r.t. larger particle size to obtain the
                        // expected sigma = 1, also converts to reduced units and puts
                        // particles back in box (whether or not PBC are used in the 
                        // following analysis)
                        boxx /= (2.0f * rTemp);
                        boxy /= (2.0f * rTemp);
                        boxz /= (2.0f * rTemp);
                        for (int i = 0; i < N; i++){
                                particles[i].x /= (2.0f * rTemp);
                                particles[i].y /= (2.0f * rTemp);
                                particles[i].z /= (2.0f * rTemp);
                                particles[i].r = rTemp;
                                particles[i].index = i;
                                particles[i].x = fmod(particles[i].x + 2 * boxx, boxx);
                                particles[i].y = fmod(particles[i].y + 2 * boxy, boxy);
                                particles[i].z = fmod(particles[i].z + 2 * boxz, boxz);
                                particles[i].type = 'a';
                        }
                        break;
                case 1: ;
                        double boxxmin = 100.0f;
                        double boxxmax = 0.0f;
                        double boxymin = 100.0f;
                        double boxymax = 0.0f;
                        double boxzmin = 100.0f;
                        double boxzmax = 0.0f;
                        for (int i = 0; i < N; i++){
                                if (fscanf(initfile, "%*s %lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].z) != 3){
                                        printf("error here 2\n");
                                        exit(0);
                                }
                                //else {printf("%lf\n", particles[i].x); printf("%lf\n", particles[i].y); printf("%lf\n", particles[i].z);}
                                particles[i].type = 'a';
                                particles[i].r = 0.5f;
                                particles[i].index = i;
                                //printf("%d\n", particles[i].index);
                                // Attempt to determine boxsize
                                if (particles[i].x < boxxmin) {boxxmin = particles[i].x;}
                                if (particles[i].x > boxxmax) {boxxmax = particles[i].x;}
                                if (particles[i].y < boxymin) {boxymin = particles[i].y;}
                                if (particles[i].y > boxymax) {boxymax = particles[i].y;}
                                if (particles[i].z < boxzmin) {boxzmin = particles[i].z;}
                                if (particles[i].z > boxzmax) {boxzmax = particles[i].z;}
                        }
                        boxx = boxxmax - boxxmin;
                        boxy = boxymax - boxymin;
                        boxz = boxzmax - boxzmin;
                        for (int i = 0; i < N; i++) {
                            Particle* p = &(particles[i]);
                            p->x -= boxxmin;
                            p->y -= boxymin;
                            p->z -= boxzmin;
                        }
                        break;
        }
        // Counts the number of read snapshots
        snapshot_count++;
        // Places cursor at the current reading position
        cursor_current = ftell(initfile);
}



void writeCoords(char *filename, int fluidlike, int snap){
/* Function:    writeCoords
 * -------------------------
 * Writes the position, radius, and type data for all N
 * particles in a .sph file; standard units are used
 * NB: Writing is done using the typesetting indicated above
 *
 * *filename:  pointer to the name of the output .sph file
 * fluidlike:  choice to print out fluid-like particles with reduced size for
 *             visualization purposes
 */
    FILE *outfile = NULL;
    outfile = fopen(filename, "a");
    if (outfile != NULL) {
        switch (output_filetype) {
            case 0:
                fprintf(outfile, "&%d\n", N);
                fprintf(outfile, "%.12lf %.12lf %.12lf\n", boxx * sigma, boxy * sigma, boxz * sigma);
                switch (fluidlike){
                    case 0:
                        for (int i = 0; i < N; i++){
                            fprintf(outfile,
                                "%c %.12lf %.12lf %.12lf %.12lf\n",
                                particles[i].type,
                                particles[i].x * sigma,
                                particles[i].y * sigma,
                                particles[i].z * sigma,
                                particles[i].r * sigma
                                );
                        }
                        break;
                    case 1:
                        for (int i = 0; i < N; i++){
                            if (particles[i].type == 'a'){
                                fprintf(outfile,
                                    "%c %.12lf %.12lf %.12lf %.12lf\n",
                                    particles[i].type,
                                    particles[i].x * sigma,
                                    particles[i].y * sigma,
                                    particles[i].z * sigma,
                                    particles[i].r * sigma / 5.0f
                                    );
                            }
                            else{
                                fprintf(outfile,
                                    "%c %.12lf %.12lf %.12lf %.12lf\n",
                                    particles[i].type,
                                    particles[i].x * sigma,
                                    particles[i].y * sigma,
                                    particles[i].z * sigma,
                                    particles[i].r * sigma
                                    );
                            }
                        }
                        break;
                }
                break;
            case 1:
                fprintf(outfile, "%d\n", N);
                fprintf(outfile, "Snap %d Box %.12lf %.12lf %.12lf\n", snap, boxx, boxy, boxz);
                for (int i = 0; i < N; i++) {
                    fprintf(outfile,
                    "%c %.12lf %.12lf %.12lf\n",
                    particles[i].type,
                    particles[i].x * sigma,
                    particles[i].y * sigma,
                    particles[i].z * sigma
                    );
                }
                break;
        }
        fclose(outfile);
    }
    else {printf("[clusterID] ERROR: Could not open output file %s\n", filename); exit(3);}
}



int buildCL(int usecase){
/* Function:   buildCL
 * --------------------
 * Builds, updates, and expands cell lists and their table according to
 * the new state of the system
 *
 * usecase:    choice of what operation to perform on CL among:
 *             + 1:    initialization of the CL table and CLs
 *                     contents
 *             + 2:    update of the CLs contents only, no
 *                     redefinition of the table, its size, or the
 *                     size of the cells
 *             + 3:    redefinition of the size of the cell, which
 *                     decreases, and update of the CLs contents; less
 *                     of the previously allocated memory is used in
 *                     this particular case
 *
 * return:      0 for normal termination
 */
        switch (usecase)
        {
                case 1:
                // sCell goes back to ~sigma
                // nCell goes up, i.e. more cells
                // more memory is allocated to accomodate for increase in number of cells
                {
                        sCellx = boxx / ((int) (boxx / sigma));
                        nCellx = (int) (boxx / sCellx);
                        sCelly = boxy / ((int) (boxy / sigma));
                        nCelly = (int) (boxy / sCelly);
                        sCellz = boxz / ((int) (boxz / sigma));
                        nCellz = (int) (boxz / sCellz);
                        nCell = nCellx * nCelly * nCellz;
                        free(CLTable);
                        CLTable = malloc(nCell * sizeof(**CLTable));
                        if (CLTable == NULL) {printf("error CL\n"); exit(EXIT_FAILURE);}
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 2: 
                // sigma < sCell < 2*sigma
                // same number of cells
                // no more memory allocated
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        break;
                }
                case 3: 
                // sCell goes back to ~sigma
                // nCell goes down, i.e. less cells
                // no more memory allocated, we juste use less of it
                {
                        for (int i = 0; i < nCell; i++)
                                CLTable[i] = NULL;
                        sCellx = boxx / ((int) (boxx / sigma));
                        nCellx = (int) (boxx / sCellx);
                        sCelly = boxy / ((int) (boxy / sigma));
                        nCelly = (int) (boxy / sCelly);
                        sCellz = boxz / ((int) (boxz / sigma));
                        nCellz = (int) (boxz / sCellz);
                        nCell = nCellx * nCelly * nCellz;
                        break;
                }
                default:
                        return 1;
        }
        for (int i = 0; i < N; i++)
        {
                Particle *one = particles+i;
                int u = one->x / sCellx;
                int v = one->y / sCelly;
                int w = one->z / sCellz;
                int index = retrieveIndex(u, v, w);
                one->CLindex = index;
                one->prev = NULL;
                one->next = CLTable[index];
                if (CLTable[index] != NULL)
                        CLTable[index]->prev = one;
                CLTable[index] = one;
        }
        return 0;
}



int retrieveIndex(int u, int v, int w){
/* Function:   retrieveIndex
 * --------------------------
 * Retrieves index of the cell (from the table of CLs) to which a
 * particle belongs based on the coordinates of the CL in the 3D
 * coordinate system of the CL table; corresponds of a flattening
 * operation of the CL index from a 3D view to standard 1D view
 *
 * u:          x-coordinate of the cell
 * v:          y-coordinate of the cell
 * w:          z-coordinate of the cell
 *
 * return:     index of the CL in the table of CLs
 */
        int index = 0;
        index = (w + 2 * nCellz) % nCellz * nCelly * nCellx
                + (v + 2 * nCelly) % nCelly * nCellx
                + (u + 2 * nCellx) % nCellx;
        return index;
}



int findClusters(void){
/* Function:   findClusters
 * ------------------------
 * Wrapper for the identification of clusters from the local BOO
 * analysis
 * Updates the table accounting for the number of clusters of a given
 * size
 *
 * return:     size of the largest cluster in the current state of the
 *             system
 */
        int maxclus;     
        build_nblist(NNidMethod);
        order = calc_order();
        connections = calc_conn(order);
        maxclus = calc_clusters(connections, order);
        return maxclus;
}



void build_nblist(int method){
/* Function:    build_nblist
 * -------------------------
 * build the list of nearest neighbours (NN) for all particles in the system
 * based on the choice of identification method
 *
 * method:      choice of NN identification method among:
 *              + 0:    SANN
 *              + 1:    12 NN
 *              + 2:    cutoff radius rc
 */
        for (int p = N-1; p > -1; p--){
                // local variables:
                int id;
                bndT* bnd;
                double d,
                       dxy,
                       dx,
                       dy,
                       dz,
                       di;

                // finding to which cell belongs the current particle
                int a = particles[p].x / sCellx;
                int b = particles[p].y / sCelly;
                int c = particles[p].z / sCellz;
        
                // initializes number of nearest neighbours for particle `p` to 0
                blist[p].n = 0;

                // defining another buffer particle
                int cellIndex;
                Particle *current = NULL;

                // saving info about neigbouring particles
                posnb neighbours[N];
                int numposnb = 0;

                // cycle over all particles in all neighboring CLs
                for (int i = a-2; i < a+3; i++){
                        for (int j = b-2; j < b+3; j++){
                                for (int k = c-2; k < c+3; k++){
                                        if ((!PBC) && ((i > nCellx) || (i < 0) || (j > nCelly) || (j < 0) || (k > nCellz) || (k < 0))) continue;
                                        cellIndex = retrieveIndex(i,j,k);
                                        current = CLTable[cellIndex];
                                        while (current != NULL && numposnb <= N){
                                                int currentID = current->index;
                                                if (p == currentID){
                                                // escapes if both particles are the same
                                                        current = current->next;
                                                        continue;
                                                }
                                                // computing distances
                                                dx = particles[p].x - current->x;
                                                dy = particles[p].y - current->y;
                                                dz = particles[p].z - current->z;
                                                // applying NIC
                                                dx = dx - boxx * rint(dx / boxx);
                                                dy = dy - boxy * rint(dy / boxy);
                                                dz = dz - boxz * rint(dz / boxz);
                                                // computes squared distance between two particles
                                                d = dx * dx + dy * dy + dz * dz;
                                                // saves neighbour data
                                                neighbours[numposnb].part = current;
                                                neighbours[numposnb].dist = sqrt(d); 
                                                neighbours[numposnb].dx = dx;                   
                                                neighbours[numposnb].dy = dy;
                                                neighbours[numposnb].dz = dz;
                                                // continues with the next neighbouring particle
                                                numposnb ++;
                                                current = current->next;
                                        }
                                }
                        }
                } 
                // sorts neighbours by increasing distanc:
                qsort(neighbours, numposnb, sizeof(posnb), compare);
                // identify NN based on method of choice 
                int m = 0;
                switch (method){
                        case 0:
                                m = 3; 
                                int done = 0;
                                while (!done){
                                        double rim = 0;
                                        for (int i = 0; i < m; i++)
                                        rim += neighbours[i].dist / (m - 2);
                                        if (rim > neighbours[m].dist) {m++;}
                                        else {done = 1;}
                                        if (m > numposnb){ 
                                                printf("[clusterID] NN algorithm did not converge!\n");
                                                m--;
                                                done = 1;
                                        }
                                }
                                for (int i = 0; i < m; i++){
                                        posnb* nb = neighbours+i;
                                        id = nb->part->index;
                                        dx = nb->dx;
                                        dy = nb->dy;
                                        dz = nb->dz;
                                        di = 1.0 / (nb->dist);
                                        bnd = &(blist[p].bnd[blist[p].n]);
                                        blist[p].n++;
                                        bnd->n = id;
                                        bnd->nz = dz * di;
                                        dxy = 1.0 / sqrt(dx * dx + dy * dy);
                                        bnd->si = dy * dxy;
                                        bnd->co = dx * dxy;
                                }
                                break;
                        case 1:
                                m = 12;
                                for (int i = 0; i < m; i++){
                                        posnb* nb = neighbours+i;
                                        id = nb->part->index;
                                        dx = nb->dx;
                                        dy = nb->dy;
                                        dz = nb->dz;
                                        di = 1.0 / (nb->dist);
                                        bnd = &(blist[p].bnd[blist[p].n]);
                                        blist[p].n++;
                                        bnd->n = id;
                                        bnd->nz = dz * di;
                                        dxy = 1.0 / sqrt(dx * dx + dy * dy);
                                        bnd->si = dy * dxy;
                                        bnd->co = dx * dxy;
                                }
                                break;
                        case 2:
                                //for (int i = 0; i < MAX_NEIGHBORS; i++){
                                for (int i = 0; i < fmin(numposnb, MAX_NEIGHBORS); i++){
                                        posnb* nb = neighbours+i;
                                        //printf("%.4lf\n", nb->dist);
                                        if (nb->dist < bndLength){
                                                //printf("added?\n");
                                                bnd = &(blist[p].bnd[blist[p].n]);
                                                blist[p].n++;
                                                id = nb->part->index;
                                                bnd->n = id;
                                                if (id > p){
                                                        //printf("added!\n");
                                                        dx = nb->dx;
                                                        dy = nb->dy;
                                                        dz = nb->dz;
                                                        di = 1.0 / (nb->dist);
                                                        bnd->nz = dz * di;
                                                        dxy = 1.0 / sqrt(dx * dx + dy * dy);
                                                        bnd->si = dy * dxy;
                                                        bnd->co = dx * dxy;
                                                }
                                        }
                                }
                                break;
                }
        }
}



int compare(const void * a, const void * b){
/* Function:    compare
 * --------------------
 * For use with qsort() function
 *
 * return:      1      if a > b
 *              0      if a = b
 *             -1      if a < b
 */
        double d = ((posnb*) a)->dist - ((posnb*) b)->dist;
        return ((0 < d) - (d < 0));
}



compl* calc_order(void){
/* Function:    calc_order
 * -----------------------
 * Computes the normalized qs (ten Wolde's local orientational order parameter) 
 *
 * return:     table of size 13*N of complex type of the normalized qs
 */
        compl* q1;
        compl* q2;
        const int l = 6;
        double temp;
        memset(orderp, (int) 0.0, sizeof(compl) * N * (l * 2 + 1));
        
        for(int i = 0; i < N; i++){
                q1 = (orderp + i * (2 * l + 1) + l);
                for(int j = 0; j < blist[i].n; j++){
                        if(blist[i].bnd[j].n > i) {
                                q2 = (orderp + blist[i].bnd[j].n * (2 * l + 1) + l);
                                compute_order(l, &(blist[i].bnd[j]), q1, q2);
                        }  
                }
        }  
        // normalize qs
        for(int i = 0; i < N; i++){
                temp = sqrt(dotprod(orderp + i * (2 * l + 1), orderp + i * (2 * l + 1), l));
                temp = 1.0 / temp;
                for(int m = -l ; m <= l; m++){
                        (*(orderp+i * (2 * l + 1) + m + l)).re *= temp;
                        (*(orderp+i * (2 * l + 1) + m + l)).im *= temp;
                }
        }
        return orderp;
}



double dotprod(compl *vec1, compl *vec2, int l){
/* Function:    dotprod
 * --------------------
 * Computes the dot product of two complex vectors
 *
 * *vec1:      pointer to a first complex vector
 * *vec2:      pointer to a second complex vector
 *
 * return:     dot product of two complex vectors
 */
        double res = 0.0f;
        for(int m = -l; m <= l; m++)
                res += (*(vec1 + m + l)).re * (*(vec2 + m + l)).re
                       + (*(vec1 + m + l)).im * (*(vec2 + m + l)).im;
        return res;
}



double sqr(double x){
/* Function:    sqr
 *  ---------------
 *  Computes the square of a double
 *
 *  x:          double to compute the square of
 *
 *  return:     x squared
 */
        return x*x;
}



void compute_order(int l, bndT *bnd, compl *res1, compl *res2){
/* Function:    order
 * ------------------
 * Computes the spherical harmonics for two neighbouring particles and feeds
 * the sum of said spherical harmonics of all neighbouring particles for a
 * given particle in calculating the local orientational order parameter
 *
 * l:           = 6
 * bnd:         table of information about the bond
 * res1:        sum of the spherical harmonics for particle 1
 * res2:        sum of the spherical harmonics for particle 2
 */
        double fc,
               p,
               f,
               s,
               r,
               sp,
               spp,
               c,
               cp,
               cpp;
        double z;
        int m = 0;
        z = bnd->nz;

        // Computes the spherical harmonics for m = 0
        p = plgndr(l,0,z);
        fc = facs(l,0);
        f = sqrt((2*l+1) * INVFPI * fc);
        r = p*f;
        (res1+0)->re += r;
        (res1+0)->im += 0;
        (res2+0)->re += r * minpow(l); // minpow(6)=1
        (res2+0)->im += 0;
        
        s=0;
        sp=0;
        c=1;
        cp=0;

        for(m = 1; m <= l; m++){
                // For m > 0
                p = plgndr(l,m,z);
                fc = facs(l,m);
                f = sqrt((2 * l + 1) * INVFPI * fc);
                r = p * f;
                // Chebyshev recursive method for computing cosine of multiple angles 
                cpp = cp;
                cp = c;
                if(m == 1){c = bnd->co;}
                else{c = 2.0 * bnd->co * cp - cpp;}
                // Chebyshev recursive method for computing sine of multiple angles 
                spp = sp;
                sp = s;
                if(m == 1){s = bnd->si;}
                else{s = 2.0 * bnd->co * sp - spp;}
                
                (res1+m)->re += r*c;
                (res1+m)->im += r*s;
                (res2+m)->re += r*c;
                (res2+m)->im += r*s;
                
                // For m < 0
                r *= minpow(m);
                (res1-m)->re += r*c;
                (res1-m)->im += -r*s;
                (res2-m)->re += r*c;
                (res2-m)->im += -r*s;
        }
}



float plgndr(int l, int m, double x){
/* Function:    plgndr
 * -------------------
 * Calculates the Legendre function P_{l,m}(x) = (1-x**2)**{m/2} (\frac{d}{dx})**m P_l(x)
 * where P_l(x) is the Legendre polynomial defined for x on [-1;1]
 *
 * l:           parameter
 * m:           parameter
 * x:           cos(theta), must be in [-1;1]
 *
 * return:      P_{l,m}(x)
 */
        // variables
        double fact,
               pll = 0.0,
               pmm,
               pmmp1,
               somx2;
        int i,
            ll;

        // checks for normal computation of Legendre polynoms
        if (m < 0 || m > l || fabs(x) > 1.0)
                printf("Bad arguments in routine plgndr %i %i %f\n", l, m, fabs(x));
        
        pmm = 1.0;
        
        if (m > 0){
                somx2 = sqrt((1.0 - x) * (1.0 + x));
                fact = 1.0;
                for (i = 1; i <= m; i++){
                        pmm *= -fact * somx2;
                        fact += 2.0;
                }
        }

        if (l == m)
                return pmm;
        else{ 
                pmmp1 = x * (2 * m + 1) * pmm;
                if (l == (m + 1))
                        return pmmp1;
                else{
                        for (ll = m + 2; ll <= l; ll++){
                                pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
                                pmm = pmmp1;
                                pmmp1 = pll;
                        }
                        return pll;
                }
        }
}



double facs(int l, int m){
/* Function:    facs
 * -----------------
 * Computes (l-m)!/(l+m)!
 *
 * l:           parameter
 * m:           parameter
 *
 * return:      (l-m)!/(l+m)!
 */
        static double* fac_table = NULL;
        int a, b;
        if(fac_table == NULL){
                fac_table = malloc((2*l+1) * sizeof(*fac_table));
                for(a = 0; a < 2*l+1; a++){
                        b = a - l;
                        fac_table[a]= exp(gammln(l - b + 1) - gammln(l + b + 1));
                }
        }
        return fac_table[m+l];
}



double gammln(double xx){
/* Function:    gammln
 * -------------------
 * Uses the Gamma function to compute factorials
 *
 * xx:          input number
 *
 * return:      factorial of xx-1
 */
        double x,
               y,
               tmp,
               ser;
        static double cof[6] = {76.18009172947146,
                                -86.50532032941677,
                                24.01409824083091,
                                -1.231739572450155,
                                0.1208650973866179e-2,
                                -0.5395239384953e-5
                                };
        int j;
        y = x = xx;
        tmp = x + 5.5;
        tmp -= (x + 0.5) * log(tmp);
        ser = 1.000000000190015;
        for (j = 0; j <= 5; j++) 
                ser += cof[j] / ++y;
        return -tmp + log(2.5066282746310005 * ser / x);
}



double minpow(int m){
/* Function:    minpow
 * -------------------
 * Returns 1.0f if m is even, -1.0f if m is odd
 *
 * m:           input integer
 *
 * return:      1.0f if m is even, -1.0f if m is odd
 */
        if((m & 1) == 1)
                return -1.0;
        else
                return 1.0;
}



int* calc_conn(compl* orderp){
/* Function:    calc_conn
 * ----------------------
 * Determines whether a particle is to be considered solid-like or not
 *
 * orderp:      table of normalized local oriental order parameters
 *
 * return:      table of int indicating whether a particle i is solid-like
 *              (conn[i]=1) or not (conn[i]=0)
 */
        int z;
        const int l = 6;
        for(int i = 0; i < N; i++){
                z = 0;
                for(int j = 0; j < blist[i].n; j++){
                        if(dotprod(orderp + i * (2 * l + 1), orderp + blist[i].bnd[j].n * (2 * l + 1), l) > bnd_cutoff){ 
                                z++;
                        }
                }
                if(z >= nbnd_cutoff){
                // particle is solid-like
                        conn[i]=1;
                } 
                else {
                // particle is fluid-like
                        conn[i]=0;
                }
                particles[i].clusterindex = 0;
        }
        return conn;
}



void setcluss(int pn, compl* orderp, int cn){
/* Function:    setcluss
 * ---------------------
 * Identifies cluster from potential solid-like particles in the list of NN
 */
        const int l = 6;
        cluss[pn] = cn;
        for(int jj = 0; jj < blist[pn].n; jj++){
                int tmp = blist[pn].bnd[jj].n;
                if(conn[tmp] != 0
                   && cluss[tmp] == 0
                   && dotprod(orderp + pn * (2 * l + 1), orderp + tmp * (2 * l + 1), l) > obnd_cutoff
                   ){
                //if(particle tmp is solid-like
                //   && particle tmp does not yet belong to a cluster
                //   && the dot product between q6(tmp) and q6(pn) > obnd_cutoff
                //   )
                //obnd_cutoff = 0.9 gives nice results
                //obnd_cutoff = 0 gives all touching nuclei as one big nuclei
                        cs++; //cluster size goes up
                        particles[tmp].clusterindex = cn;
                        setcluss(tmp, orderp, cn);
                }  
        }
        if (highlightAll && conn[pn] == 1) particles[pn].type = 'a' + cn%27;
}



int calc_clusters(int* conn, compl* orderp){
/* Function:   calc_clusters
 * -------------------------
 *
 * conn:        table of int indicating whether a particle i is solid-like
 *              (conn[i]=1) or not (conn[i]=0)
 * orderp:      table of size 13*N of complex type of the normalized qs
 *
 * return:      size of the largest cluster in the current state of the
 *              system
 */
        // VARIABLES
        int cn = 1;             //cluster id number
        int big = 0;            //largest cluster size
        int bc = -1;            //largest cluster id number
        int unread = 1;


        // BODY
        for(int i = 0; i < N; i++){
                cluss[i] = 0;           
                size[i] = 0;
        }  

        for(int i = 0; i < N; i++){
                cs = 0;
                if(conn[i] == 1 && cluss[i] == 0){
                //if(particle i is solid-like && particle i does not yet belong to a cluster)
                        cs++;
                        setcluss(i, orderp, cn); 
                        size[cn] = cs;
                        if(cs > big){
                        //identifies largest cluster
                                big = cs;
                                bc = cn;
                        }
                        particles[i].clusterindex = cn;
                        if (highlightAll && conn[i] == 1) particles[i].type = 'a' + cn%27;
                        cn++;
                }
        }

        if (!highlightAll){
        // highlight largest cluster only
                for (int i = 0; i < N; i++){
                        if (particles[i].clusterindex == bc)
                                particles[i].type = 'b';        
                }
        }
        
        logBook->bookmark = 0;
        for (int i = 0; i < cn; i++){
                if (size[i] != 0 && logBook->bookmark == 0){
                // first passage
                        logBook->size[logBook->bookmark] = size[i];
                        logBook->qt[logBook->bookmark] = 1;
                        logBook->bookmark++;
                        Nn[size[i]] += 1.0f / (double) N;
                }
                else if (size[i] != 0){
                        unread = 1;
                        for (int j = 0; j < logBook->bookmark; j++){
                                if (size[i] == logBook->size[j]){
                                        logBook->qt[j]++;
                                        unread=0;
                                }
                        }
                        if (unread){
                                logBook->size[logBook->bookmark] = size[i];
                                logBook->qt[logBook->bookmark] = 1;
                                logBook->bookmark++;
                        }
                        Nn[size[i]] += 1.0f / (double) N;
                }
        }
        
        return big;
}



void saveLogBook(void){
/* Function:    saveLogBook
 * ------------------------
 * Save logs of identified clusters to a file. Logs are not saved in increasing
 * cluster size.
 */
        FILE* saveclusters = fopen(clusterlog_filename, "a");
        if (saveclusters != NULL){
                fprintf(saveclusters, "&%d\n", logBook->bookmark);
                for (int i = 0; i < logBook->bookmark; i++){
                        fprintf(saveclusters, "%d\t%d\n", logBook->size[i], logBook->qt[i]);
                }
                fclose(saveclusters);
        }
}



void saveDGndata(void){
/* Function:    saveDGndata
 * ------------------------
 * Saves DGn data to a file.
 * The calculation is not rigourous for experimental data that was collected
 * with a fluctuating total number of particles between snapshots (this can be
 * due to poor particles tracking).
 */
        FILE* DGnfile = fopen(DGndata_filename, "a");
        if (DGnfile != NULL){
                for (int i = 0; i < N; i++){
                        if (Nn[i] > 0.0f){
                                fprintf(DGnfile, "%d\t%.12lf\n", i, -log(Nn[i] / (double) snapshot_count));
                        }
                }
                fclose(DGnfile);
        }
}





//    MAIN
// ==========
int main(int argc, char *argv[]) {
        // CPU TIME MEASUREMENT | START
        struct timespec beginc, endc;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginc);
        // WALL TIME MEASUREMENT | START
        struct timespec beginw, endw;
        clock_gettime(CLOCK_REALTIME, &beginw);



        // PARSING INPUT
        parse_input(argc, argv);



        // VERBOSE
        printf("[clusterID] Higlighting clusters based on the following parametrization:\n * SOURCE:\t\t%s\n * bnd_cutoff:\t\t%.2lf\n * nbnd_cutoff:\t\t%d\n * obnd_cutoff:\t\t%.2lf\n * NN id method:\t", input_filename, bnd_cutoff, nbnd_cutoff, obnd_cutoff);
        if (NNidMethod == 0) printf("SANN\n");
        else if (NNidMethod == 1) printf("12NN\n");
        else if (NNidMethod == 2) printf("cutoff radius = %.2lf\n", bndLength);
        if (savelogs) printf(" * Writing of cluster logs enabled\n");
        if (energy) printf(" * Calculation of DGn enabled\n");
        printf("[clusterID] Starting...\n");



        // OUTPUT
        if (save_movie) {
            if (output_filetype == 0) {
                if (sprintf(movie_filename, "./movie.sph") < 0){
                    printf("[clusterID] Could not write string. Exit.\n");
                    exit(3);
                }
            }
            else if (output_filetype == 1) {
                if (sprintf(movie_filename, "./movie.xyz") < 0) {
                    printf("[clusterID] Could not write string. Exit.\n");
                    exit(3);
                }
            }
            FILE* writefile = fopen(movie_filename, "w+");
            if (writefile != NULL) fclose(writefile);
        }

        if (savelogs){
                if (sprintf(clusterlog_filename, "./logs.txt") < 0){
                        printf("[clusterID] Could not write string. Exit.\n");
                        exit(0);
                }
                FILE* writefile = fopen(clusterlog_filename, "w+");
                if (writefile != NULL){fclose(writefile);}
        }

        if (energy){
                if (sprintf(DGndata_filename, "./DGn.txt") < 0){
                        printf("[clusterID] Could not write string. Exit.\n");
                        exit(0);
                }
                FILE* writefile = fopen(DGndata_filename, "w+");
                if (writefile != NULL){fclose(writefile);}
        }



        // BODY
        FILE* initfile = fopen(input_filename, "r");
        snap_count = 0;
        if (initfile != NULL){
                initReading(initfile);
                while (cursor_current < cursor_end-10){
                        readCoords(initfile);
                        buildCL(1);
                        findClusters();
                        if (save_movie) writeCoords(movie_filename, 1, snap_count);
                        if (savelogs) saveLogBook();
                        printf("[clusterID] Snapshot %d done.\n", snap_count);
                        snap_count++;
                        // Freeing allocated memory before allocating again
                        free(particles);
                        free(connections);
                        free(order);
                        free(cluss);
                        free(size);
                        for (int i = 0; i < N; i++)
                                free(blist[i].bnd);
                        free(blist);
                        free(logBook->size);
                        free(logBook->qt);
                        free(logBook);
                }
        }
        fclose(initfile);
        if (energy) saveDGndata();

       

        // END
        free(Nn);
        printf("[clusterID] Program terminated normally\n * Read snapshots:\t%d\n", snapshot_count);
        if (save_movie) printf(" * Produced\t\t%s\n", movie_filename);
        if (savelogs) printf(" * Produced\t\t%s\n", clusterlog_filename);
        if (energy) printf(" * Produced\t\t%s\n", DGndata_filename);
        //  CPU TIME MEASUREMENT | END
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endc);
        long secondsc = endc.tv_sec - beginc.tv_sec;       
        long nanosecondsc = endc.tv_nsec - beginc.tv_nsec;
        double elapsedc = secondsc + nanosecondsc * 1e-9;
        printf(" * Elapsed CPU time:\t%.3fs\n", elapsedc);
        // WALL TIME MEASUREMENT | END
        clock_gettime(CLOCK_REALTIME, &endw);
        long secondsw = endw.tv_sec - beginw.tv_sec;       
        long nanosecondsw = endw.tv_nsec - beginw.tv_nsec;
        double elapsedw = secondsw + nanosecondsw * 1e-9;
        printf(" * Elapsed wall time:\t%.3fs\n", elapsedw);

        return 0;
}
