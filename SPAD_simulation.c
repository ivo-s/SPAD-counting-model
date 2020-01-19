/*
 * SPAD_simulation.c
 * 
 * Simulation of SPAD detections - outputs counting statistics.
 * Needs an afterpulsing generation table - see comment on line 262
 */

#include<pthread.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>
#include<unistd.h>
#include<inttypes.h>

// multithreading
#define NUM_THREADS 3

// 1 MICROSECOND is the unit of time here
#define WINDOW 10.	// time window for count statistics

#define DEADTIME 0.024	// recovery time
#define AP_MEAN 0.0171	// afterpulse probability
#define TWILIGHT_COEFF 0.0032	// afterpulse-mean linear coefficient.
								// Input reset time in microseconds.

// for nonlinear twilight pulsing
//~ #define TWILIGHT_PROB_EXPLICIT 0.002

// number of detections to simulate
#define SAMPLES ((long long) 1e9)

// allocation constants
#define DEFAULT_AP_QUEUE_SIZE 10
#define AP_MAX_ORDER 10
#define HIST_MAX_SIZE 1000

// number of random generators
#define NUMBER_OF_PRNGS 4
#define PRNG_BUFSIZE 256

// normalization constant for generating random floats
#define RAND_M_F (((double)RAND_MAX)+1.)

// afterpulse table floating point format
typedef float APfloat;

// struct for passing to child threads
typedef struct {
	int threadNum;
	APfloat *APdist;	// AP distribution
	long dist_size;		// AP dist. size
	double mean;		// mean incident rate
	unsigned int seed[NUMBER_OF_PRNGS];  // seeds for random generators
	volatile unsigned long long samples; // sample monitoring
	volatile int go;	// for telling the thread to stop
	unsigned long long *returnHist;
	double simTime;
} childArg;

// struct for afterpulses waiting to happen
typedef struct {
	double *queue;
	int len;		// population length
	int alloc_len;	// allocation length
} APqueue;

void RemoveAPfromQueue(APqueue *apq, int i) {
	for (apq->len--; i < apq->len; i++) {
		apq->queue[i] = apq->queue[i+1];
	}
}

void AddAPtoQueue(APqueue *apq, double t) {
	if (apq->len == apq->alloc_len) {	// check allocated space
		void *realloc_ptr;
		realloc_ptr = realloc( apq->queue,
			(apq->alloc_len+10)*sizeof(*apq->queue) );
		if (realloc_ptr == NULL) {
			fprintf(stderr, "ERROR: cannot reallocate AP queue.\n");
			return;
		}
		apq->alloc_len += 10;
		apq->queue = realloc_ptr;
	}
	apq->queue[apq->len++] = t;
}

// the main simulation function run by each thread
void * generateHist(void *arg) {
	
	unsigned long long i, j;

	childArg *a = (childArg *) arg;
	unsigned long long hist[HIST_MAX_SIZE];
	for (i=0; i<HIST_MAX_SIZE; i++) hist[i] = 0LL;
	
#ifndef TWILIGHT_PROB_EXPLICIT
	double twilightProb = a->mean*TWILIGHT_COEFF;
#else
	double twilightProb = TWILIGHT_PROB_EXPLICIT;
#endif	
	
	// AP queue allocation
	APqueue apq;
	apq.len = 0;
	apq.alloc_len = DEFAULT_AP_QUEUE_SIZE;
	apq.queue = malloc(DEFAULT_AP_QUEUE_SIZE*sizeof(double));
	if (apq.queue == NULL) {
		fprintf(stderr, "ERROR: cannot allocate memory for AP queue\n");
	}

	// Poisson CDF for generating afterpulses
	double poisson[AP_MAX_ORDER];
	poisson[0] = exp(-AP_MEAN);
	for (i=1; i<AP_MAX_ORDER; i++) {
		poisson[i] = (poisson[i-1]*AP_MEAN)/i;
	}
	for (i=1; i<AP_MAX_ORDER; i++) {
		poisson[i] += poisson[i-1];
	}

	// declaration of random generators
	// see manpage random_r(3)
	struct random_data buf[NUMBER_OF_PRNGS];
	memset(&buf[0], 0, NUMBER_OF_PRNGS*sizeof(struct random_data));
	char statebuf[NUMBER_OF_PRNGS][PRNG_BUFSIZE];
	memset(&statebuf[0][0], 0, NUMBER_OF_PRNGS*PRNG_BUFSIZE);
	int32_t r_int;

	double t=0.;	// time
	double r, rRad;	// for storing random numbers
	int nAP;		// number of afterpulses
	int indexAP;	// for generating afterpulse delay
	int count=0;	// main detection counter

	double tOverflow = 0.;	// time overflow

	// PRNG initialisation
	for (i=0; i<NUMBER_OF_PRNGS; i++) {
		if (initstate_r(	a->seed[i], &statebuf[i][0], PRNG_BUFSIZE,
							&buf[i]	))
		{
			fprintf(stderr, "initstate_r ERROR\n");
		}
	}
	
	// MAIN LOOP
	// each loop generates the time of one detection
	for (i=0; a->go; ++i) {

		t += DEADTIME; // add after each detection

		// clear old afterpulses
		j=0;
		while (j<apq.len) {
			if (apq.queue[j] <= t) {
				RemoveAPfromQueue(&apq, j);
			} else j++;
		}
		
		// create afterpulses as a point process

		// generate nAP as a Poissonian variable
		random_r(&buf[0], &r_int);	// PRNG 1
		r = r_int/RAND_M_F;
		for (nAP=0; nAP < AP_MAX_ORDER; nAP++) {
			if (r < poisson[nAP]) break;
		}
		
		// generate delay for each afterpulse
		for (j=0; j<nAP; j++) {
			random_r(&buf[1], &r_int);	// PRNG 2
			indexAP = (int) ((r_int/RAND_M_F)*(a->dist_size));
			AddAPtoQueue(&apq, t + a->APdist[indexAP]);
		}
		
		random_r(&buf[2], &r_int);	// PRNG 3
		r = r_int/RAND_M_F;
		// only add time if there is no twilight pulse
		if (r > twilightProb) {
			random_r(&buf[3], &r_int);	// PRNG 4
			rRad = r_int/RAND_M_F;
			// initialize time to the next Poissonian event
			t = t - log(1.-rRad)/(a->mean);
			// search for afterpulses that come before
			for (j=0; j<apq.len; j++) {
				if (apq.queue[j] < t) {
					t = apq.queue[j];
				}
			}
		}
		
		// check time window and wrap around
		
		while (t >= WINDOW) {
			t -= WINDOW;
			for (j=0; j<apq.len; j++) {
				apq.queue[j] -= WINDOW;
			}
			tOverflow += WINDOW;
			// increment histogram
			if (count < HIST_MAX_SIZE) hist[count]++;
			count = 0;
		}
		count++;
		
		// update samples every now and then
		if (i % 1000000 == 0) {
			a->samples = i;
		}
	}
	
	// simulation statistics
	a->simTime = tOverflow + t;
	a->samples = i;

	// return the result
	a->returnHist = (unsigned long long *) malloc(sizeof(unsigned long long)*HIST_MAX_SIZE);
	memcpy(a->returnHist, hist, sizeof(unsigned long long)*HIST_MAX_SIZE);
	
	free(apq.queue);
	
	return 0;

}

int main(int argc, char *argv[]) {
	
	if (argc < 4) {
		printf("ARGUMENTS: <AP dist table> <incident mean rate> <output file>\n");
		return 0;
	}
	
	// threading variables
	pthread_t threads[NUM_THREADS];
	int thread_res;
	childArg thread_args[NUM_THREADS];
	unsigned long long overallProgress;
	// while thread_args can be close together in memory and are
	// accessed by different threads, the writing is kept to a minimum
	
	unsigned long long i;
	int j;
	
	// afterpulsing distribution
	APfloat *APdist;
	long dist_size;
	
	// histogram variables
	int lastElement;
	unsigned long long hist[HIST_MAX_SIZE];
	for (i=0; i<HIST_MAX_SIZE; i++) hist[i] = 0LL;
	
	unsigned long int counts = 0;
	double simTime = 0.;

	// import afterpulse time distribution table

	// needs to be a sequence of APfloats distributed so that if
	// we pick an element uniformly at random, the element would be
	// distributed in time according to the AP distribution
	
	// in other words, an inverse CDF of the AP delay
	
	FILE *d = fopen(argv[1], "rb");
	if (d == NULL) {
		fprintf(stderr, "Error opening file %s\n", argv[1]);
		return 0;
	}
	fseek(d, 0, SEEK_END);
	APdist = (APfloat *) malloc(ftell(d));
	if (APdist == NULL) {
		fprintf(stderr, "Failed to allocate memory\n");
		fclose(d);
		return 0;
	}
	if (ftell(d) % sizeof(APfloat) != 0) {
		fprintf(stderr, "File does not contain an integer number of floats\n");
		fclose(d);free(APdist);
		return 0;
	}
	dist_size = ftell(d)/sizeof(APfloat);
	fseek(d, 0, SEEK_SET);
	if (fread(APdist, sizeof(APfloat), dist_size, d) == 0) {
		fprintf(stderr, "Error reading from %s\n", argv[1]);
		fclose(d);free(APdist);
		return 0;
	}
	fclose(d);
	
	// open output file
	FILE *f = fopen(argv[3],"w");
	if (f == NULL) {
		fprintf(stderr, "Cannot open output file %s\n", argv[3]);
		free(APdist);
		return 0;
	}

	// random seed
	srandom(time(NULL));
	
	// initialize thread variables
	for (i=0; i<NUM_THREADS; i++) {
		thread_args[i].threadNum = i;
		thread_args[i].APdist = APdist;
		thread_args[i].dist_size = dist_size;
		sscanf(argv[2], "%lf", &thread_args[i].mean); // read mean rate
		for (j=0; j<NUMBER_OF_PRNGS; j++) {	// initialize random seeds
			thread_args[i].seed[j] = random();
		}
		thread_args[i].samples = 0LL;
		thread_args[i].go = 1;
	}
	
	// start threads
	for (i=0; i<NUM_THREADS; i++) {
		thread_res = pthread_create(&threads[i], NULL, &generateHist, &thread_args[i]);
		if (thread_res != 0) {
			fprintf(stderr, "Starting thread %llu resulted in %d\n", i, thread_res);
		}
	}
	
	// waiting loop
	while (1) {
		overallProgress = 0LL;
		// sum the numbers of samples
		for (i=0; i<NUM_THREADS; i++) {
			overallProgress += thread_args[i].samples;
		}
		printf("%2.2lf %%\r", (float) overallProgress/SAMPLES*100);
		fflush(stdout);
		// if enough samples have been generated, break
		if (overallProgress >= SAMPLES) break;
		usleep(100000);
	}
	
	// set the red light for the loops to finish
	for (i=0; i<NUM_THREADS; i++) {
		thread_args[i].go = 0;
	}
	
	for (i=0; i<NUM_THREADS; i++) {
		thread_res = pthread_join(threads[i], NULL);
		if (thread_res != 0) {
			fprintf(stderr, "Joining thread %llu resulted in %d\n", i, thread_res);
		}
		// add the simulated histogram to the total
		for (j=0; j<HIST_MAX_SIZE; j++) {
			hist[j] += thread_args[i].returnHist[j];
		}
		free(thread_args[i].returnHist);
		// simulation time and counts
		simTime += thread_args[i].simTime;
		counts += thread_args[i].samples;
	}
	
	// find the index of the last non-zero element
	for (lastElement=HIST_MAX_SIZE-1; lastElement>0; lastElement--) {
		if (hist[lastElement] > 0) break;
	}

	printf("Counts:    %lu\nTime:      %lf\nDet. rate: %lf\n", counts, simTime, counts/simTime);
	
#ifndef TWILIGHT_PROB_EXPLICIT
	double twilightProb = thread_args[0].mean*TWILIGHT_COEFF;
#else
	double twilightProb = TWILIGHT_PROB_EXPLICIT;
#endif
	
	// print file header with parameters
	fprintf( f, "WINDOW %lf\nRATE_PER_TIMEUNIT %.16lf\nDEADTIME %lf\n"
			"AFTERPULSE_MEAN %lf\nTWILIGHT_PROB %lf\n", WINDOW,
			thread_args[0].mean, DEADTIME, AP_MEAN, twilightProb);

	// print the histogram
	for (i=0; i<=lastElement; i++) {
		if (i>0) putc(' ', f);
		fprintf(f, "%llu", hist[i]);
	}
	
	fclose(f);
	free(APdist);

	return 0;

}
