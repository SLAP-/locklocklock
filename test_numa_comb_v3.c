//The queue lock version of measuring the local, remote delays and handover probabilities
#define _GNU_SOURCE
#include <stdio.h>
#include <numa.h>
#include <pthread.h>
#include <stdlib.h>
#include <sys/types.h> //for pthread_spinlock_t
#include <unistd.h>
#include <errno.h>
#include "tsc.h"
#include "cs_queue.c"

#define N_THREADS 8
#define CPU_FREQ 2270000000
#define EXPERIMENT_TIME_IN_SEC 15
#define MAX_N_ACCESS 100000
#define S_TO_N 1000000000

//function prototypes
static int cmp_timestamp(const void *ptr1, const void *ptr2);
void set_affinity(pthread_t thread, const cpu_set_t *cpuset);
void *work(void *thread_arg);
double abs_value(double a);
int is_on_same_node(int i, int j, int n, int left, int right);

//#define DEBUG
double time_in_cs[N_THREADS];
int num_access_each_thread[N_THREADS];
volatile int start_work_flag = 0;
pthread_t threads[N_THREADS];
cpu_set_t cpuset[N_THREADS];
int access_count = 0;

typedef struct
{
	int thread_id;
	struct random_data* rand_states;
	int arrival_lambda;
	pthread_spinlock_t *spinlock_ptr;
}thread_params;

typedef struct
{
	uint64_t ts;
	int id;
}timestamp;

timestamp g_tss[MAX_N_ACCESS*N_THREADS]; //global timestamps

static int cmp_timestamp(const void *ptr1, const void *ptr2)
{
	timestamp *a = (timestamp *)ptr1;
	timestamp *b = (timestamp *)ptr2;
	if(a->ts > b->ts)
		return 1;
	else if(a->ts == b->ts)
		return 0;
	else
		return -1;
}
void set_affinity(pthread_t thread, const cpu_set_t *cpuset)
{
	int s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), cpuset);
	if (s != 0)
	{
		perror("set affinity error \n");
		exit(-1);
	}
}
void *work(void *thread_arg)
{
	while(!start_work_flag)
		;
	thread_params *thread_param = (thread_params*)thread_arg;
	int tid = thread_param->thread_id;
	//int arrival_lambda = thread_param->arrival_lambda;
	double sevice_time = 0.001;
	//int initial_loop_flag = 1;
	timestamp tss[MAX_N_ACCESS]; //timsstampes
	int count = 0;
	set_affinity(threads[tid], &cpuset[tid]);
	uint64_t start_time = read_tsc();
	
	while(read_tsc() - start_time < EXPERIMENT_TIME_IN_SEC*CPU_FREQ) 
	{
		uint64_t start = read_tsc();
		while(read_tsc() - start < sevice_time*CPU_FREQ/2)
			;
		uint64_t getlock_time = read_tsc();
		#ifdef DEBUG
		printf("Thread %d trying to get the lock at %lu \n",tid,getlock_time);
		#endif
		impl_enter_critical(tid);
		//pthread_spin_lock(thread_param->spinlock_ptr);	
		#ifdef DEBUG
		printf("Thread %d got the lock at %lu \n",tid,read_tsc());
		#endif
		/************The critical section***************/
		uint64_t get_in_cs_time = read_tsc();
		//get_lock_time[access_count] = get_in_cs_time;
		tss[count++].ts = get_in_cs_time;
		while(read_tsc() - get_in_cs_time < sevice_time*CPU_FREQ)
		{
		#ifdef DEBUG
			if(if_print_cs_info == 1) 
			{
				printf("Thread %d in cs.\n",tid);
				if_print_cs_info = 0;
			}
		#endif
			;
		}
		//cs_order[access_count] = tid; //Record my id in the cs order array
		//printf("access count: %d\n",access_count);
		//release_lock_time[access_count] = read_tsc();
		tss[count++].ts =  read_tsc();
		//access_count++;
		/************End the critical section***************/
		impl_exit_critical(tid);
		//pthread_spin_unlock(thread_param->spinlock_ptr);
		uint64_t end_cs_time = read_tsc();
#ifdef DEBUG
		printf("Thread %d released the lock at %lu \n",tid,end_cs_time);
#endif
		num_access_each_thread[tid]++;
		time_in_cs[tid] +=  (end_cs_time - getlock_time)/(double)CPU_FREQ;
	}
	//when the experiment is done, write to the global timestamp array
	pthread_spin_lock(thread_param->spinlock_ptr);	
		int i;
		for(i = 0; i < count; i++)
		{
			g_tss[access_count].ts = tss[i].ts;
			g_tss[access_count].id = tid;
			access_count++ ;
		}
	pthread_spin_unlock(thread_param->spinlock_ptr);
	return NULL;
}

double abs_value(double a)
{
	if(a > 0)
		return a;
	else
		return -a;
}
int is_on_same_node(int i, int j, int n, int left, int right)
{
	if(right == 0)
		return 1;
	if((i < left && j >= left) || (i >= left && j < left))
		return 0;
	else
		return 1;
}
int main(int argc, char *argv[])
{
	int node_id = 0;
	int arrival_lambda = 10;
	int thread_cpu_map[N_THREADS];
	int i,j,k;
	int n_threads;
	int n_left;
	int n_right;
	int next_index_left = 3;
	int next_index_right = 7;
	float local_square = 0.0, remote_square = 0.0;
	/*
	   pthread_spinlock_t *spinlock_ptr = malloc(sizeof(pthread_spinlock_t));
	   if(spinlock_ptr == NULL) //error handling of the malloc of the spinlock
	   {
	   printf("malloc of spinlock failed.\n");
	   }
	   else
	   {
	   printf("malloc of spinlock succeeded.\n");
	   }
	   free(spinlock_ptr);
	 */
	/*****************make sure the number of arguments is correct and get the n_threads, n_left and n_right*/
	if(argc < 4)
	{
		printf("Usage: ./test_numa_comb number_of_threads number_of_threads_on_node0 number_of_threads_on_node1\n");
		exit(-1);
	}
	n_threads = atoi(argv[1]);
	n_left = atoi(argv[2]);
	n_right = atoi(argv[3]);
	/*******************Set the thread_cpu_map according to the n_left and n_right*/
	printf("n_threads: %d, n_left: %d, n_right: %d\n",n_threads,n_left,n_right);
	for(i = 0; i < n_left; i++)
	{
		thread_cpu_map[i] = next_index_left;
		next_index_left--;
	}
	for(i = n_left; i < n_threads; i++)
	{
		thread_cpu_map[i] = next_index_right;
		next_index_right--;
	}
	for(i = 0; i < n_threads; i++)
	{
		printf("Thread %d is on cpu %d\n",i,thread_cpu_map[i]);
	}

	thread_params para[n_threads]; //The parameters to pass to the threads

	//printf("The return value of numa_get_run_node_mask(void) is %d\n",numa_get_run_node_mask());
	//printf("The return value of numa_max_node(void) is %d\n",numa_max_node());
	//numa_tonode_memory((void *)spinlock_ptr,sizeof(pthread_spinlock_t),node_id); //This doesn't work

	//initilize the spinlock pointer and put it on a specific node
	pthread_spinlock_t *spinlock_ptr = numa_alloc_onnode(sizeof(pthread_spinlock_t),node_id);
	if(spinlock_ptr == NULL) //error handling of the allocating of a spinlock pointer on a specific node
	{
		printf("alloc of spinlock on a node failed.\n");
	}
	else
	{
		printf("alloc of spinlock on a node succeeded.\n");
	}
	//initlize spinlock
	pthread_spin_init(spinlock_ptr,0);
	impl_init(n_threads);
	//create the threads
	for(i = 0; i < n_threads; i++)
	{
		int rc;
		para[i].thread_id = i;
		para[i].arrival_lambda = arrival_lambda;
		para[i].spinlock_ptr = spinlock_ptr;
		CPU_ZERO(&cpuset[i]);
		CPU_SET(thread_cpu_map[i],&cpuset[i]);
		rc = pthread_create(&threads[i],NULL,work,(void*)&para[i]);
		if(rc)
		{
			printf("ERROR: return code from pthread_create() is %d for thread %d \n",rc,i);
			exit(-1);
		}
		start_work_flag = 1; 
	}
	for(i = 0; i < n_threads; i++)
	{
		pthread_join(threads[i],NULL);
	}
	/*
	for(i = 0; i < n_threads; i++)
	{
		printf("The time to get one lock for thread %d is : %.9f\n",i,time_in_cs[i]/num_access_each_thread[i]);
		printf("The number of lock accesses for thread %d is : %d\n",i,num_access_each_thread[i]);
	}
	*/
	//double diff = abs_value(time_in_cs[0]/num_access_each_thread[0] - time_in_cs[1]/num_access_each_thread[1]);
	//printf("The difference of the time to get one lock is : %.9f (%f cycles)\n",diff,diff*CPU_FREQ); //this is assuming there are only two processors (needs to be changed if there are more)

	//printf("printing the get lock and release lock time:\n");
	qsort((void*)g_tss,(size_t)access_count,(size_t)sizeof(timestamp),cmp_timestamp);
	/*
	for (i = 0; i < access_count; i++)
		printf("%lu with id %d\n",g_tss[i].ts,g_tss[i].id);
	*/
	/*
	for (i = 1; i < access_count - 2; i+=2)
	{
		double d = (double)(g_tss[i+1].ts - g_tss[i].ts);
		//printf("%d release to %d require: %.9f seconds (%.9f nanoseconds,%d cycles)\n",g_tss[i].id,g_tss[i+1].id,d/CPU_FREQ,(d*S_TO_N)/CPU_FREQ,d);
	}
	*/
	int cs_order[access_count/2];
	for(i = 0; i < access_count/2; i++)
	{
		cs_order[i] = g_tss[i*2].id;
		//printf("%d in cs\n",cs_order[i]);
	}
	int cs_matrix[n_threads][n_threads];
	uint64_t delay_matrix[n_threads][n_threads];
	float prob_matrix[n_threads][n_threads];
	int local_count2 = 0, remote_count2 = 0;
	for(i = 0; i < n_threads; i++)
	{
		for(j = 0; j < n_threads ; j++)
		{
			cs_matrix[i][j] = 0;
			delay_matrix[i][j] = 0;
			prob_matrix[i][j] = 0.0;

		}
	}
	for(i = 0; i < n_threads; i++)
	{
		for(j = 0; j < n_threads; j++)
		{
			for(k = 0; k < access_count/2 -1 ; k++)
			{
				if(cs_order[k] == i && cs_order[k+1] == j)
				{
					cs_matrix[i][j]++;
					delay_matrix[i][j] += (g_tss[2*k+2].ts - g_tss[2*k+1].ts); 
					if(is_on_same_node(i, j, n_threads, n_left, n_right))
					{
						printf("local_delay: %lu\n",(g_tss[2*k+2].ts - g_tss[2*k+1].ts));
						local_square += (g_tss[2*k+2].ts - g_tss[2*k+1].ts)*(g_tss[2*k+2].ts - g_tss[2*k+1].ts); 
						local_count2++;
					}
					else
					{
						printf("remote_delay: %lu\n",(g_tss[2*k+2].ts - g_tss[2*k+1].ts));
						remote_square += (g_tss[2*k+2].ts - g_tss[2*k+1].ts)*(g_tss[2*k+2].ts - g_tss[2*k+1].ts); 
						remote_count2++;
					}
				}
			}
		}
	}
	int num_access[n_threads];
	for(i = 0; i < access_count/2 -1; i++)
	{
		for(j = 0; j < n_threads; j++)
		{

			if(cs_order[i] ==  j)
			{
				num_access[j]++;
			}
		}
	}
	/*
	for(i = 0; i < n_threads; i++)
		printf("num_access[%d]:%d\n",i,num_access[i]);
	*/

	for(i = 0; i < n_threads; i++)
	{
		for(j = 0; j < n_threads ; j++)
		{
			prob_matrix[i][j] = (float)cs_matrix[i][j]/(float)num_access[i];
			//printf("%d followed by %d: %f,",i,j, prob_matrix[i][j]);
		}
		printf("\n");
	}
	printf("Delay in cycles:\n");
	for(i = 0; i < n_threads; i++)
	{
		for(j = 0; j < n_threads ; j++)
		{
			printf("delay of %d to %d: %f,",i,j,delay_matrix[i][j]/(float)cs_matrix[i][j]);
		}
		printf("\n");
	}
	//print the intra-core and inter-core delay
	//thread 0 - n_left -1 are on the left core, n_left to n_threads are on the right core
	uint64_t local_delay = 0, remote_delay = 0;
	int local_count = 0, remote_count = 0;
	float local_prob = 0.0, remote_prob = 0.0;

	for(i = 0; i < n_threads; i++)
	{
		for(j = 0; j < n_threads; j++)
		{
			if (j == i)
				continue;
			if(is_on_same_node(i, j, n_threads, n_left, n_right))
			{
				printf("%d and %d on the same node\n",i,j);
				local_delay += delay_matrix[i][j];
				local_count += cs_matrix[i][j];
				local_prob += prob_matrix[j][i];
			}
			else
			{
				printf("%d and %d not the same node\n",i,j);
				remote_delay += delay_matrix[i][j];
				remote_count += cs_matrix[i][j];
				remote_prob += prob_matrix[j][i];
			}
		}
	}

	float local = (float)local_delay/(local_count);
	float remote = (float)remote_delay/(remote_count);
	printf("local delay: %f, remote_delay: %f, local_count: %d, remote_count: %d\n",(float)local_delay/(local_count),(float)remote_delay/(remote_count),local_count,remote_count);
	printf("local prob:%f, remote prob: %f\n",local_prob/n_threads, remote_prob/n_threads);
	printf("local delay variance:%f, remote delay variance: %f\n",local_square/local_count - local*local, remote_square/remote_count - remote*remote);
	printf("local count2: %d, remote_count2:%d\n",local_count2, remote_count2);
	pthread_spin_destroy(spinlock_ptr);
	numa_free((void *)spinlock_ptr,sizeof(pthread_spinlock_t));
	pthread_exit(NULL);
	return 0;
}
