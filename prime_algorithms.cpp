#include <stdio.h>	//printf
#include <time.h>
#include <stdint.h>	//uint.._t
#include <math.h>	//log()
#include <iostream>
#include <vector>

#define pierre_dussart_min 355991
#define miller_rabin_max 4759123141

using namespace std;

void print_stats(string algorithm, uint64_t counter, clock_t time){
	algorithm += ':';
	cout << algorithm << string(16 - algorithm.length(), ' ')<< counter;
	printf(" primes in %lf\n", ((double)time)/CLOCKS_PER_SEC);
}

uint64_t pierre_dussart(uint64_t x){
	double lnx = log(x);
	return((x/lnx)*(1 + 1/lnx + 2.51/(lnx * lnx)));
}

uint64_t prime_estimator(uint64_t i){
	if(i >= pierre_dussart_min)
		return(pierre_dussart(i));
	return (i/2 + 1);
}

void deterministic(uint64_t current, uint64_t range){
	if(current > range || range >= -2)
		return;

	uint64_t counter = 0;
	clock_t time = clock();

	if(range >= 2){
		uint32_t i;
		if(current <= 2){
			current = 2;
			counter++;
			current++;
		}
		for(;current <= 3 && current <= range; current++)
			for(i = 2; current % i != 0; i++)
				if(i*i > current){
					counter++;
					break;
				}
		if(current % 2 == 0 && current <= range)
				current++;
		for(;current <= range; current += 2)
			for(i = 3; current % i != 0; i += 2)
				if(i*i > current){
					counter++;
					break;
				}
	}
	time = clock() - time;
	print_stats("Deterministic", counter, time);
}

void sieve_of_atkin(uint64_t minimum, uint64_t range){
	if(minimum > range || range >= -2)
		return;

	uint64_t current = 0, counter = 0, last = 0;
	vector<uint64_t> prime_vector(prime_estimator(range), 0);
	prime_vector[0] = 2;

	clock_t time = clock();
	if(range >= 2){
		uint32_t i;
		if(current <= 2){
			current = 2;
			if(minimum <= 2)
				counter++;
			if(current <= range)
				current++;
		}
		for(;current <= 3 && current <= range; current++)
			for(i = 0;current % prime_vector[i] != 0; i++)
				if(current < prime_vector[i] * prime_vector[i]){
					prime_vector[last+1] = current;
					last ++;
					counter +=  (current >= minimum);
					break;
				}
		for(;current < minimum;current++)
			for(i = 0;current % prime_vector[i] != 0; i++)
				if(current < prime_vector[i] * prime_vector[i]){
					prime_vector[last+1] = current;
					last ++;
					break;
				}
		if(current % 2 == 0 && current <= range)
			current++;
		for(; current <= range; current += 2)
			for(i = 0; current % prime_vector[i] != 0; i++)
				if(current < (prime_vector[i] * prime_vector[i])){
					prime_vector[last+1] = current;
					last ++;
					counter ++;
					break;
				}
	}
	time = clock() - time;
	print_stats("Sieve of Atkin", counter, time);
}

uint64_t modulo(uint64_t p, uint64_t k, uint64_t m){
	if(k == 0)
		return 1;
	else if(k & 1)
		return((p * modulo(p,k-1,m)) % m);
	else
		return(modulo((p*p) % m, k >> 1, m));
}

uint64_t ml_calc(uint8_t a, uint64_t r,uint64_t d,uint64_t current){
	uint64_t x = modulo(a,d,current);
	if(x == 1||x == current - 1)
		return 1;
	for(uint64_t j = 0;j < r-1 && r > 1;j++){
		x = modulo(x, 2, current);
		if (x == current - 1)
			return 1;
	}
	return 0;
}

void miller_rabin(uint64_t current, uint64_t range){
	if(current > range || range > miller_rabin_max)
		return;

	uint64_t counter = 0, r, d;
	clock_t time = clock();
	if(range >= 2){
		if(current <= 2){
			current = 2;
			counter++;
			if(current <= range)
				current++;
		}
		for(; current <= 3 && current <= range; current++)
			for(uint32_t i = 2; current % i != 0; i++)
				if(i*i > current){
					counter++;
					break;
				}
		if(current % 2 == 0 && current <= range)
			current++;
		for(; current <= range; current += 2){
			if(current == 7 || current == 61){
				counter++;
				continue;
			}
			for(r = 0;((current - 1) >> r) % 2 == 0;r++);
			d = ((current - 1) >> r);
			counter += (ml_calc(2, r, d, current) && ml_calc(7, r, d, current) && ml_calc(61, r, d, current));
		}
	}
	time = clock() - time;
	print_stats("Miller-Rabin", counter, time);
}

int main(int argc, char *argv[]){
		cout << endl;
	if(argc != 3)
		cout << "Parameters: min max" << endl;
	else{
		uint64_t min = atoi(argv[1]);
		uint64_t max = atoi(argv[2]);
		cout << "Checking range [" << min << ", " << max << "]" << endl;
		cout << "Prime estimation: " << prime_estimator(max) << endl;
		miller_rabin(min, max);
		sieve_of_atkin(min, max);
		deterministic(min, max);
	}
	cout << endl;
	return 0;
}