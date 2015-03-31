/*
 * Copyright(C) 2014 Pedro H. Penna <pedrohenriquepenna@gmail.com>
 * 
 * friendly-numbers.c - Friendly numbers kernel.
 */

#include <global.h>
#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <util.h>
#include <timer.h>
#include "fn.h"

/*
 * Computes the Greatest Common Divisor of two numbers.
 */
static int gcd(int a, int b)
{
  int c;
  
  /* Compute greatest common divisor. */
  while (a != 0)
  {
     c = a;
     a = b%a;
     b = c;
  }
  
  return (b);
}

/*
 * Some of divisors.
 */
static int sumdiv(int n)
{
	int sum;    /* Sum of divisors. */
	int factor; /* Working factor.  */
	
	sum = 1 + n;
	
	/* Compute sum of divisors. */
	for (factor = 2; factor < n; factor++)
	{
		/* Divisor found. */
		if ((n%factor) == 0)
			sum += factor;
	}
	
	return (sum);
}

/*
 * Computes friendly numbers.
 */
int friendly_numbers(int start, int end) 
{
	int n;        /* Divisor.                    */
	int *num;     /* Numerator.                  */
	int *den;     /* Denominator.                */
	int range;    /* Range of numbers.           */
	int i, j;     /* Loop indexes.               */
	int nfriends; /* Number of friendly numbers. */
	int tid;
	uint64_t t1, t2, *total;
	
	nfriends = 0;
	range = end - start + 1;
	
	num = smalloc(sizeof(int)*range);
	den = smalloc(sizeof(int)*range);
	
	total = scalloc(nthreads, sizeof(uint64_t));
	
	/* Compute abundances. */
	#pragma omp parallel private(i, j, tid, t1, t2, n)
	{
		tid = omp_get_thread_num();
		
		#pragma omp for schedule(static)
		for (i = start; i <= end; i++) 
		{
			t1 = timer_get();
			
			j = i - start;
							
			num[j] = sumdiv(i);
			den[j] = i;
				
			n = gcd(num[j], den[j]);
			num[j] /= n;
			den[j] /= n;
			
			t2 = timer_get();
			
			total[tid] += timer_diff(t1, t2);
		}
		
		/* Check friendly numbers. */
		#pragma omp for reduction(+:nfriends) schedule(static)
		for (i = 1; i < range; i++)
		{
			t1 = timer_get();
			
			for (j = 0; j < i; j++)
			{
				/* Friends. */
				if ((num[i] == num[j]) && (den[i] == den[j]))
					nfriends++;
			}
			
			t2 = timer_get();
			total[tid] += timer_diff(t1, t2);
		}
	}
	
	for (i = 0; i < nthreads; i++)
		fprintf(stderr, "  thread %d: %f\n", i, total[i]*MICROSEC);

	free(num);
	free(den);
	
	return (nfriends);
}
